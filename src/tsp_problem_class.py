"""TSP-specific copy of the DDRO branch-and-cut problem.

Connectivity cuts are mandatory feasibility constraints.  They are separated
at every integer lazy-callback candidate and deliberately do not share the
intersection-cut per-node limit or root-only setting.
"""

import time

import bnc_problem_class as base
import globals as g
from parse_sp import parse_sp_instance
from python_to_cplex_c_api import cplex_c_api_wrapper as cpx
from tsp_topology import (
    TSPEdge,
    TSPTopology,
    canonical_component_shores,
    connected_components,
    is_hamiltonian_tour,
)


TSP_RESULT_FIELDS = base.RESULT_FIELDS + (
    "connectivity_cut_num",
    "connectivity_time",
    "connectivity_callback_count",
    "tour_valid",
)


class ConnectivitySeparator:
    """Separate global cutset inequalities at an integer TSP candidate."""

    def __init__(self, model):
        self.num_vertices = model.tsp_topology.num_vertices
        self.edges = model.tsp_topology.edges
        self.edge_col_indices = model.edge_col_indices
        self.seen_cuts = set()

    def components_from_solution(self, solution):
        values = [solution[index] for index in self.edge_col_indices]
        return connected_components(self.num_vertices, self.edges, values)

    def separate(self, callback):
        started = time.time()
        base.tracker.connectivity_callback_count += 1

        node_lp = cpx.CPXgetcallbacknodelp(
            g.env, callback.cbdata, callback.wherefrom
        )
        num_cols = cpx.CPXgetnumcols(g.env, node_lp)
        # Read the actual integer candidate supplied to the lazy callback.
        # CPXgetx(node_lp) can instead expose the current node LP solution when
        # the candidate originated from a heuristic.
        solution = cpx.CPXgetcallbacknodex(
            g.env, callback.cbdata, callback.wherefrom, 0, num_cols - 1
        )
        components = self.components_from_solution(solution)

        if len(components) == 1:
            base.tracker.connectivity_time += time.time() - started
            return 0

        added = 0
        for shore in canonical_component_shores(self.num_vertices, components):
            if shore in self.seen_cuts:
                continue

            shore_set = set(shore)
            cut_indices = []
            for edge, col_index in zip(self.edges, self.edge_col_indices):
                if (edge.first in shore_set) != (edge.second in shore_set):
                    cut_indices.append(col_index)

            if not cut_indices:
                raise RuntimeError(
                    f"No graph edge crosses disconnected TSP shore {shore}."
                )

            cpx.CPXcutcallbackadd(
                g.env,
                callback.cbdata,
                callback.wherefrom,
                len(cut_indices),
                2.0,
                "G",
                cut_indices,
                [1.0] * len(cut_indices),
                0,
            )
            self.seen_cuts.add(shore)
            added += 1

        base.tracker.connectivity_cut_count += added
        base.tracker.connectivity_time += time.time() - started
        return added


class TSPConnectivityCallback(cpx.CutCallback):
    """Mandatory connectivity-only lazy callback used by BB-CONN."""

    def __init__(self, model):
        super().__init__()
        self.connectivity = ConnectivitySeparator(model)

    def __call__(self):
        base.tracker.callback_count += 1
        added = self.connectivity.separate(self)
        return cpx.CPX_CALLBACK_SET if added else cpx.CPX_CALLBACK_DEFAULT


class TSPIntersectionCutCallback(base.IntersectionCutCallback):
    """Connectivity first, followed by DDRO intersection separation."""

    def __init__(self, model):
        super().__init__(model)
        self.connectivity = ConnectivitySeparator(model)

    def __call__(self):
        added = self.connectivity.separate(self)
        if added:
            # This candidate is not a tour.  Cut it off before spending time on
            # the DDRO separation MIP.  Connectivity has no max-cuts limit.
            base.tracker.callback_count += 1
            return cpx.CPX_CALLBACK_SET
        return super().__call__()


class TSPBnCProblem(base.BnCProblem):
    """BOBILib-compatible BnC problem with mandatory TSP connectivity."""

    def __init__(self, con, tra):
        self.tsp_config = con
        super().__init__(con, tra)

    def generate_problem_data(self):
        super().generate_problem_data()
        if self.tsp_config.instance_type != "bobilib":
            raise ValueError("The TSP BnC copy currently requires MPS/AUX BOBILib input.")
        if not self.tsp_config.tsp_source_file:
            raise ValueError(
                "The conversion manifest must identify the original SP source file."
            )

        sp_instance = parse_sp_instance(self.tsp_config.tsp_source_file)
        self.tsp_topology = TSPTopology(
            sp_instance.num_vertices,
            tuple(
                TSPEdge(f"x_{edge.index}", edge.first, edge.second)
                for edge in sp_instance.edges
            ),
        )
        self.tsp_topology.validate()
        if not self.ul_var_names:
            raise ValueError("The BOBILib parser did not preserve upper-level names.")

        upper_name_to_col = {
            name: index for index, name in enumerate(self.ul_var_names)
        }
        missing = [
            edge.variable
            for edge in self.tsp_topology.edges
            if edge.variable not in upper_name_to_col
        ]
        if missing:
            preview = ", ".join(missing[:8])
            suffix = " ..." if len(missing) > 8 else ""
            raise ValueError(
                f"TSP edge variables are not upper-level MPS variables: "
                f"{preview}{suffix}"
            )

        self.edge_col_indices = [
            upper_name_to_col[edge.variable] for edge in self.tsp_topology.edges
        ]
        for edge, index in zip(self.tsp_topology.edges, self.edge_col_indices):
            if (
                self.lb[index] < -self.tsp_config.tolerance
                or self.ub[index] > 1.0 + self.tsp_config.tolerance
            ):
                raise ValueError(
                    f"TSP edge variable {edge.variable} must have bounds within [0, 1]."
                )

    def set_callbacks(self):
        if base.config.intersection_cuts:
            self.lazy_cb = TSPIntersectionCutCallback(self)
            self.incumbentcb = base.MyIncumbentCallback(self)
            cpx.CPXsetlazyconstraintcallbackfunc(g.env, self.lazy_cb)
            cpx.CPXsetincumbentcallbackfunc(g.env, self.incumbentcb)

            if base.config.separation == "fractional":
                # Fractional connectivity separation would require min-cut and
                # is intentionally outside this integer lazy-callback version.
                self.user_cb = base.IntersectionCutCallback(self)
                cpx.CPXsetusercutcallbackfunc(g.env, self.user_cb)
        elif base.config.branchandbound:
            self.lazy_cb = TSPConnectivityCallback(self)
            self.incumbentcb = base.MyIncumbentCallback(self)
            cpx.CPXsetlazyconstraintcallbackfunc(g.env, self.lazy_cb)
            cpx.CPXsetincumbentcallbackfunc(g.env, self.incumbentcb)
        else:
            raise ValueError(
                "TSP copy supports only 'branchandbound' (BB-CONN) and "
                "'intersection' (BC-INT-CONN)."
            )

    def finalize_result(self, result):
        """Attach TSP statistics and independently validate any incumbent tour."""
        result["problem"] = "TSP"
        result["method"] = (
            "BC-INT-CONN" if base.config.intersection_cuts else "BB-CONN"
        )
        result["connectivity_cut_num"] = base.tracker.connectivity_cut_count
        result["connectivity_time"] = base.tracker.connectivity_time
        result["connectivity_callback_count"] = (
            base.tracker.connectivity_callback_count
        )

        tour_valid = ""
        if result["has_incumbent"]:
            solution = cpx.CPXgetx(g.env, g.m, 0, self.no_of_vars - 1)
            edge_values = [solution[index] for index in self.edge_col_indices]
            tour_valid = int(is_hamiltonian_tour(
                self.tsp_topology.num_vertices,
                self.tsp_topology.edges,
                edge_values,
            ))
            if not tour_valid:
                raise RuntimeError(
                    "CPLEX returned an incumbent that is not a connected "
                    "degree-two Hamiltonian tour. Check the MPS degree constraints "
                    "and edge-variable mapping."
                )
        result["tour_valid"] = tour_valid
        return result
