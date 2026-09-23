##################################################################
# This file is part of the code used for the computational study #
# in the paper                                                   #
#                                                                #
#  "Branch-and-Cut for Mixed-Integer Linear                      # 
#   Decision-Dependent Robust Optimization"                      #
#                                                                #
# by Henri Lefebvre, Martin Schmidt, Simon Stevens,              #
# and Johannes Thürauf (2026).                                   #
##################################################################

# Global imports
import os
from cplex.exceptions import CplexError
import numpy as np
import time
import sys

# Local imports
import globals as g
from ddro_parser import *

# Import CPLEX C API wrapper
from python_to_cplex_c_api import cplex_c_api_wrapper as cpx


RESULT_FIELDS = (
    "problem",
    "instance",
    "n",
    "density",
    "replicate",
    "formulation",
    "linearcons",
    "method",
    "projected",
    "status",
    "solver_status",
    "has_incumbent",
    "incumbent",
    "bound",
    "gap",
    "gap_percent",
    "solve_time",
    "bound_time",
    "model_time_limit",
    "total_time",
    "nodes_e",
    "cut_num",
    "cut_time",
    "callback_count",
    "cone_coefficient_time",
    "lower_level_constraints",
    "objective_scale",
    "raw_incumbent",
    "raw_bound",
)


class BnCProblem:
    """
    Base class for solving DDRO problems with BnC.
    """

    def __init__(self, con, tra):
        # Initialize config and tracker
        global config
        config = con

        global tracker
        tracker = tra

        # Initialize verbose print functionality
        global v_print
        v_print = config._v_print

        # Create an empty model to start with
        self.create_empty_model()

        # Generate problem data
        self.generate_problem_data()

        # Fill the model with the data from the instance
        if config.projected:
            v_print(1, "Filling projected model...")
            self.fill_model_projected()
        else:
            v_print(1, "Filling unprojected model...")
            self.fill_model_unprojected()

        # Write the model to an LP file for inspection
        if config.write_lps:
            cpx.CPXwriteprob(g.env, g.m, "initial_model.lp")
            cpx.CPXwriteprob(g.env, g.m, "initial_model.mps")

    #######################################################
    # Class Methods
    #######################################################

    def generate_problem_data(self):
        """Generate problem data structures."""
        if config.instance_type == "knapsack":
            parser = KNAPSACKParser(config.instance_file, config)
            problem_data = parser.parse_knapsack_instance()
        elif config.instance_type == "bobilib":
            parser = MPSAUXParser(config.instance_file, config)
            problem_data = parser.parse_bobilib_instance()

        # Upper level data (certain)
        self.c = problem_data["c"]
        self.h = problem_data["h"]
        self.A = problem_data["A"]
        self.gamma = problem_data["gamma"]

        # Upper level data (uncertain)
        self.beta = problem_data["beta"]
        self.a = problem_data["a"]
        self.b = problem_data["b"]
        self.bilinearities = problem_data["bilinearities"]

        # Lower level data
        self.D = problem_data["D"]
        self.C = problem_data["C"]
        self.alpha = problem_data["alpha"]
        
        # Bounds 
        self.lb = problem_data["lb"]
        self.ub = problem_data["ub"]

        # Problem dimensions
        self.I = problem_data["I"]  # indices for variables parameterizing uncertainty set
        self.J = problem_data["J"]  # indices for variables affected by uncertainty
        self.K = problem_data["K"]  # indices for non-affected variables
        
        self.n = problem_data["n"]  # number of upper-level variables
        self.p = problem_data["p"]  # number of variables parameterizing uncertainty set
        self.q = problem_data["q"]  # number of variables affected by uncertainty
        self.r = problem_data["r"]  # number of non-affected variables

        # Metadata
        self.n_lower_level_constraints = problem_data["n_lower_level_constraints"]
        self.instance_file = problem_data["instance_file"]

    def create_empty_model(self):
        """Create an empty CPLEX model."""
        g.env = cpx.CPXopenCPLEX()
        g.m = cpx.CPXcreateprob(g.env, "BnC Problem Model")

    def fill_model_projected(self):
        """Fill the model with the projected formulation data."""
        # Store number of variables
        self.numcols = self.n + 1
        if config.interdiction_cuts:
            self.numcols += self.p + 1
        
        # Add decision variables (x and y variables)
        var_names = [f"xI_{i}" for i in range(self.p)]  + [f"xJ_{i}" for i in range(self.q)] + [f"y_{i}" for i in range(self.r)]
        var_types = "I" * len(var_names)

        # Add eta variable
        var_names.append("eta")
        var_types += "C"

        # Add auxiliary variables for mccormick
        if config.interdiction_cuts:
            for i in range(self.p):
                var_names += [f"z_{i}"]
                var_types += "B"

        # Set objective coefficients
        obj_coeffs = [0.0] * (self.n)
        for i in range(self.p):
            obj_coeffs[self.I[i]] = self.h[i]  # h_i for x_i
        for i in range(self.q):
            obj_coeffs[self.J[i]] = self.c[i]  # c_i for y_i
        for i in range(self.r):
            obj_coeffs[self.K[i]] = self.c[i]  # Coefficient for non-affected variables y_i
        obj_coeffs.append(0.0)  # Coefficient for eta

        # Set bounds
        lb = []
        ub = []
        for i in range(self.p + self.q + self.r):  # bounds for x and y
            lb.append(self.lb[i])
            ub.append(self.ub[i])
        lb_eta, ub_eta = self.get_eta_bounds() # Upper and lower bound for eta
        lb.append(lb_eta)  # Lower bound for eta
        ub.append(ub_eta)  # Upper bound for eta

        # Set coefficients and bounds for auxiliary interdiction variables
        if config.interdiction_cuts:
            for i in range(self.n + 1, self.n + 1 + self.p):  # bounds for z
                obj_coeffs.append(0.0)
                lb.append(0.0)
                ub.append(1.0)
        
        # Add variables to the model
        cpx.CPXnewcols(
            g.env, g.m, len(obj_coeffs), obj_coeffs, lb, ub, var_types, var_names
        )
        # Store eta variable index for easy access
        self.eta_idx = self.n
        self.no_of_vars = len(var_names)

        # Change objective sense
        cpx.CPXchgobjsen(g.env, g.m, cpx.CPX_MIN)
        
        # Add constraints
        self.add_upper_level_constraints()
        self.add_epigraph_constraint()

    def get_eta_bounds(self):
        """Compute a valid lower and upper bound for eta variable."""
        if config.instance_type == "knapsack":
            return 0.0, cpx.CPX_INFBOUND
        
        lb_env = cpx.CPXopenCPLEX()
        lb_m = cpx.CPXcreateprob(lb_env, "LB Model")
        
        
        num_vars = self.r + self.p + self.q
        # Add decision variables (x and y variables)
        var_names = [f"x_{i}" for i in range(self.p)] + [f"y_{i}" for i in range(self.r)] + [f"y_{i}" for i in range(self.q)]
        var_types = "I" * len(var_names)
        
        
        # Set objective coefficients
        obj_coeffs = [0.0] * (num_vars)
        for i in range(self.p):
            obj_coeffs[i] = 0.0  # Coefficient for x_i  
        for i in range(self.r):
            obj_coeffs[self.p + i] = -self.a[i]  # c_i for y_i
        for i in range(self.q):
            obj_coeffs[self.p + self.r + i] = -self.h[i]  # c_i for y_i
        
        # Set bounds
        lb = self.lb
        ub = self.ub

        # Add variables to the model
        cpx.CPXnewcols(
            lb_env, lb_m, num_vars, obj_coeffs, lb, ub, var_types, var_names
        )
        
        # Change objective sense
        cpx.CPXchgobjsen(lb_env, lb_m, cpx.CPX_MIN)
        
        # Add upper level constraints 
        for i in range(len(self.A)):
            constraint_indices = []
            constraint_coeffs = []
            constraint_senses = []
            for j in range(len(self.A[i])):
                if self.A[i][j] != 0:
                    constraint_indices.append(j)
                    constraint_coeffs.append(self.A[i][j])
            constraint_senses = ["G"]
            # Add the constraint to the model
            cpx.CPXaddrows(
                lb_env,
                lb_m,
                0,
                1,
                len(constraint_coeffs),
                [self.gamma[i]],
                constraint_senses,
                [0],
                constraint_indices,
                constraint_coeffs,
                None,
                [f"upper_level_constr_{i}"],
            )
        if config.write_lps:
            cpx.CPXwriteprob(lb_env, lb_m, "lb_model.lp")
        
        # Optimize
        cpx.CPXmipopt(lb_env, lb_m)
        lb_eta_one = cpx.CPXgetobjval(lb_env, lb_m)
        
        cpx.CPXchgobjsen(lb_env, lb_m, cpx.CPX_MAX)
        # Optimize
        cpx.CPXmipopt(lb_env, lb_m)
        lb_eta_two = cpx.CPXgetobjval(lb_env, lb_m)
        return min(lb_eta_one, -lb_eta_one, lb_eta_two, -lb_eta_two), max(lb_eta_one, -lb_eta_one, lb_eta_two, -lb_eta_two)

    def fill_model_unprojected(self):
        """Fill the model with the unprojected formulation data."""
        # Number of variables
        self.numcols = self.n + 2*(self.r + self.q)

        # Add decision variables (x and y variables)
        var_names = (
            [f"xI_{i}" for i in range(self.p)]
            + [f"xJ_{i}" for i in range(self.q)]
            + [f"y_{i}" for i in range(self.r)]
            + [f"u_{i}" for i in range(self.r + self.q)]
            + [f"r_{i}" for i in range(self.r + self.q)]
        )
        var_types = "I" * (len(var_names))
        
        # Add eta variable
        var_names.append("eta")
        var_types += "C" 
        self.eta_idx = len(var_names) - 1

        # Add auxiliary variables for interdiction mccormick
        if config.interdiction_cuts:
            for i in range(self.p):
                var_names += [f"z_{i}"]
                var_types += "B"
            self.numcols += self.p

        # Calculate objective coefficients
        obj_coeffs = [0.0] * (self.numcols)
        for i in range(self.p):
            obj_coeffs[self.I[i]] = self.h[i]  # h_i for x_i
        for i in range(self.q):
            obj_coeffs[self.J[i]] = self.c[i]  # c_i for y_i
        for i in range(self.r):
            obj_coeffs[self.K[i]] = self.c[i]  # Coefficient for non-affected variables y_i
        obj_coeffs.append(0.0)  # Coefficient for eta
        
        # Set bounds
        lb = []
        ub = []
        for i in range(self.p + self.q + self.r):  # bounds for xI, xJ, y
            lb.append(self.lb[i])
            ub.append(self.ub[i])
        for i in range(self.q + self.r):  # bounds for u
            lb.append(self.lb[self.p + i])
            ub.append(self.ub[self.p + i])
        for i in range(self.q + self.r):  # bounds for r
            lb.append(0.0)
            ub.append(1.0)
        lb_eta, ub_eta = self.get_eta_bounds() # Upper and lower bound for eta
        lb.append(lb_eta)  # Lower bound for eta
        ub.append(ub_eta)  # Upper bound for eta

        # Set coefficients and bounds for auxiliary interdiction variables
        if config.interdiction_cuts:
            for i in range(2*self.p):  # bounds for z
                obj_coeffs.append(0.0)
                lb.append(0.0)
                ub.append(1.0)

        self.no_of_vars = len(var_names)
        
        # Add variables to the model
        cpx.CPXnewcols(
            g.env, g.m, len(var_names), obj_coeffs, lb, ub, var_types, var_names
        )
        
        # Change objective sense
        cpx.CPXchgobjsen(g.env, g.m, cpx.CPX_MIN)
        
        v_print(1, f"Filled unprojected model with {self.numcols} variables.")
        # Add constraints
        self.add_upper_level_constraints()
        self.add_lower_level_constraints()
        self.add_epigraph_constraint()

        # Add McCormick constraints for unprojected formulation
        self.add_unprojected_mccormick_constraints()

    def add_unprojected_mccormick_constraints(self):
        """Add McCormick constraints for bilinear terms r_k = x_k * u_k."""
        # For each k from 0 to min(p,q)-1, add McCormick constraints for r_k = x_k * u_k
        for k in range(min(self.p, self.q)):
            x_idx = self.J[k]  # x_k variable index
            u_idx = self.n + k  # u_k variable index
            r_idx = self.n + self.p + k  # r_k variable index

            # Constraint 1: r_k <= x_k
            constraint_indices = [r_idx, x_idx]
            constraint_coeffs = [1.0, -1.0]
            constraint_senses = ["L"]
            cpx.CPXaddrows(
                g.env,
                g.m,
                0,
                1,
                len(constraint_coeffs),
                [0.0],
                constraint_senses,
                [0],
                constraint_indices,
                constraint_coeffs,
                None,
                [f"mccormick1_r{k}"],
            )

            # Constraint 2: r_k <= u_k
            constraint_indices = [r_idx, u_idx]
            constraint_coeffs = [1.0, -1.0]
            constraint_senses = ["L"]
            cpx.CPXaddrows(
                g.env,
                g.m,
                0,
                1,
                len(constraint_coeffs),
                [0.0],
                constraint_senses,
                [0],
                constraint_indices,
                constraint_coeffs,
                None,
                [f"mccormick2_r{k}"],
            )

            # Constraint 3: r_k >= x_k + u_k - 1
            # Rearranged: -r_k + x_k + u_k <= 1
            constraint_indices = [r_idx, x_idx, u_idx]
            constraint_coeffs = [-1.0, 1.0, 1.0]
            constraint_senses = ["L"]
            cpx.CPXaddrows(
                g.env,
                g.m,
                0,
                1,
                len(constraint_coeffs),
                [1.0],
                constraint_senses,
                [0],
                constraint_indices,
                constraint_coeffs,
                None,
                [f"mccormick3_r{k}"],
            )

    def add_upper_level_constraints(self):
        """Add upper-level constraints to the model."""
        for i in range(len(self.A)):
            constraint_indices = []
            constraint_coeffs = []
            constraint_senses = []
            for j in range(len(self.A[i])):
                if self.A[i][j] != 0:
                    constraint_indices.append(j)
                    constraint_coeffs.append(self.A[i][j])
            constraint_senses = ["G"]
            # Add the constraint to the model
            cpx.CPXaddrows(
                g.env,
                g.m,
                0,
                1,
                len(constraint_coeffs),
                [self.gamma[i]],
                constraint_senses,
                [0],
                constraint_indices,
                constraint_coeffs,
                None,
                [f"upper_level_constr_{i}"],
            )

    def add_epigraph_constraint(self):
        """Add capacity constraints for the nominal knapsack model."""
        # Add capacity constraint
        constraint_indices = []
        constraint_coeffs = []
        constraint_senses = []

        # Capacity constraint: sum(a_i*x_i) + eta <= beta
        if self.bilinearities:
            for i in range(self.q):
                constraint_indices.append(self.J[i])
                constraint_coeffs.append(self.a[i])
            constraint_indices.append(self.eta_idx)  # Add eta variable index
            constraint_coeffs.append(1.0)  # Coefficient for eta
            constraint_senses.append("L")
        else:
            # TODO switch sign of constraint_coeffs depending on min or max in lower level
            for i in range(self.r):
                constraint_indices.append(self.p+i)
                constraint_coeffs.append(self.a[i])
            constraint_indices.append(self.eta_idx)  # Add eta variable index
            constraint_coeffs.append(-1.0)  # Coefficient for eta
            constraint_senses.append("L")

        # Add the constraint to the model
        cpx.CPXaddrows(
            g.env,
            g.m,
            0,
            1,
            len(constraint_coeffs),
            [self.beta],
            constraint_senses,
            [0],
            constraint_indices,
            constraint_coeffs,
            None,
            ["unc_constr_proj"],
        )

        # Add McCormick constraints for interdiction cuts
        if config.interdiction_cuts:
            self.add_mccormick_constraints()
    
    def add_mccormick_constraints(self):
        """Add McCormick constraints for bilinear terms z_k = x_k * y_k."""
        # For each k from 0 to min(p,q)-1, add McCormick constraints for z_k = x_k * y_k
        for k in range(min(self.p, self.q)):
            x_idx = self.I[k]  # x_k variable index
            y_idx = self.J[k]  # y_k variable index
            z_idx = self.n + 1 + k  # z_k variable index
            if config.projected == False:
                z_idx += self.q + self.p

            # Constraint 1: z_k <= x_k
            constraint_indices = [z_idx, x_idx]
            constraint_coeffs = [1.0, -1.0]
            constraint_senses = ["L"]
            cpx.CPXaddrows(
                g.env,
                g.m,
                0,
                1,
                len(constraint_coeffs),
                [0.0],
                constraint_senses,
                [0],
                constraint_indices,
                constraint_coeffs,
                None,
                [f"mccormick1_z{k}"],
            )

            # Constraint 2: z_k <= y_k
            constraint_indices = [z_idx, y_idx]
            constraint_coeffs = [1.0, -1.0]
            constraint_senses = ["L"]
            cpx.CPXaddrows(
                g.env,
                g.m,
                0,
                1,
                len(constraint_coeffs),
                [0.0],
                constraint_senses,
                [0],
                constraint_indices,
                constraint_coeffs,
                None,
                [f"mccormick2_z{k}"],
            )

            # Constraint 3: z_k >= x_k + y_k - 1
            # Rearranged: -z_k + x_k + y_k <= 1
            constraint_indices = [z_idx, x_idx, y_idx]
            constraint_coeffs = [-1.0, 1.0, 1.0]
            constraint_senses = ["L"]
            cpx.CPXaddrows(
                g.env,
                g.m,
                0,
                1,
                len(constraint_coeffs),
                [1.0],
                constraint_senses,
                [0],
                constraint_indices,
                constraint_coeffs,
                None,
                [f"mccormick3_z{k}"],
            )

    def add_lower_level_constraints(self):
        """Add lower-level constraints to the model."""
        # Add lower-level capacity constraint
        constraint_indices = []
        constraint_coeffs = []
        constraint_senses = []
        for j in range(len(self.C)):
            constraint_indices = []
            constraint_coeffs = []
            for i in range(len(self.C[j])):
                if self.C[j][i] != 0:
                    constraint_indices.append(i)
                    constraint_coeffs.append(self.C[j][i])
            for i in range(len(self.D[j])):
                if self.D[j][i] != 0:
                    constraint_indices.append(self.q + self.p + self.r + i)
                    constraint_coeffs.append(self.D[j][i])
            constraint_senses.append("G")
            # Add the constraint to the model
            cpx.CPXaddrows(
                g.env,
                g.m,
                0,
                1,
                len(constraint_coeffs),
                [self.alpha[j]],
                constraint_senses,
                [0],
                constraint_indices,
                constraint_coeffs,
                None,
                [f"ll_constraint_{j}"],
            )

    def set_callbacks(self):
        """Set callbacks based on configuration."""
        if config.intersection_cuts:
            self.lazy_cb = IntersectionCutCallback(self)
            self.user_cb = IntersectionCutCallback(self)
            self.incumbentcb = MyIncumbentCallback(self)

            if config.separation == "fractional":
                cpx.CPXsetusercutcallbackfunc(g.env, self.user_cb)
            cpx.CPXsetlazyconstraintcallbackfunc(g.env, self.lazy_cb)
            cpx.CPXsetincumbentcallbackfunc(g.env, self.incumbentcb)

        elif config.interdiction_cuts:
            self.lazy_cb = InterdictionCutCallback(self)
            cpx.CPXsetlazyconstraintcallbackfunc(g.env, self.lazy_cb)

        elif config.no_good_cuts:
            self.lazy_cb = NoGoodCutCallback(self)
            self.user_cb = NoGoodCutCallback(self)
            self.incumbentcb = MyIncumbentCallback(self)

            if config.separation == "fractional":
                cpx.CPXsetusercutcallbackfunc(g.env, self.user_cb)
            cpx.CPXsetlazyconstraintcallbackfunc(g.env, self.lazy_cb)
            cpx.CPXsetincumbentcallbackfunc(g.env, self.incumbentcb)

        elif config.branchandbound:
            self.incumbentcb = MyIncumbentCallback(self)
            cpx.CPXsetincumbentcallbackfunc(g.env, self.incumbentcb)

    def set_priorities(self):
        """Set variable priorities to guide branching decisions."""
        cnt = self.p + self.q + 1  # Number of variables in the problem
        indices = range(self.p + self.q + 1)
        priority = [2] * self.p + [1] * self.q + [0]  # Priorities for x_i, y_i, and eta
        direction = [-1] * len(priority)  # Directions for the variables

        # Set priorities for the variables
        cpx.CPXcopyorder(g.env, g.m, cnt, indices, priority, direction)

    def print_solution(self):
        """Print the available solution and return a CSV-ready result row."""
        solver_status = cpx.CPXgetstat(g.env, g.m)
        solution_count = cpx.CPXgetsolnpoolnumsolns(g.env, g.m)
        has_incumbent = solution_count > 0

        status_names = {
            cpx.CPXMIP_OPTIMAL: "OPTIMAL",
            cpx.CPXMIP_OPTIMAL_TOL: "OPTIMAL_TOL",
            cpx.CPXMIP_INFEASIBLE: "INFEASIBLE",
            cpx.CPXMIP_INForUNBD: "INFEASIBLE_OR_UNBOUNDED",
            cpx.CPXMIP_UNBOUNDED: "UNBOUNDED",
            cpx.CPXMIP_TIME_LIM_FEAS: "TIME_LIMIT_FEASIBLE",
            cpx.CPXMIP_TIME_LIM_INFEAS: "TIME_LIMIT_NO_INCUMBENT",
            cpx.CPXMIP_MEM_LIM_FEAS: "MEMORY_LIMIT_FEASIBLE",
            cpx.CPXMIP_MEM_LIM_INFEAS: "MEMORY_LIMIT_NO_INCUMBENT",
            cpx.CPXMIP_NODE_LIM_FEAS: "NODE_LIMIT_FEASIBLE",
            cpx.CPXMIP_NODE_LIM_INFEAS: "NODE_LIMIT_NO_INCUMBENT",
            cpx.CPXMIP_ABORT_FEAS: "ABORTED_FEASIBLE",
            cpx.CPXMIP_ABORT_INFEAS: "ABORTED_NO_INCUMBENT",
        }
        status = status_names.get(solver_status, f"CPLEX_STATUS_{solver_status}")

        objective = float("nan")
        gap = float("nan")
        if has_incumbent:
            objective = cpx.CPXgetobjval(g.env, g.m)
            gap = cpx.CPXgetmiprelgap(g.env, g.m)
            print(f"Objective value: {objective}")

            if config.verbose > 0:
                solution = cpx.CPXgetx(g.env, g.m, 0, self.no_of_vars - 1)
                print("Solution:")
                for i in range(self.p):
                    print(f"  xI_{i} = {solution[self.I[i]]}")
                for i in range(self.q):
                    print(f"  xJ_{i} = {solution[self.J[i]]}")
                for i in range(self.r):
                    r_idx = self.p + self.q + i
                    print(f"  y_{i} = {solution[r_idx]}")
                print(f"  eta = {solution[self.eta_idx]}")

                if config.projected == False:
                    for i in range(self.q):
                        u_idx = self.n + i
                        print(f"  u_{i} = {solution[u_idx]}")
                    for i in range(self.q):
                        r_idx = self.n + self.p + i
                        print(f"  r_{i} = {solution[r_idx]}")

                if config.interdiction_cuts:
                    for k in range(min(self.p, self.q)):
                        z_idx = self.n + 1 + k
                        print(f"  z_{k} = {solution[z_idx]}")
        else:
            print("No incumbent solution was found.")

        try:
            best_bound = cpx.CPXgetbestobjval(g.env, g.m)
        except RuntimeError:
            best_bound = float("nan")

        # Print additional statistics
        if config.interdiction_cuts:
            method = "interdiction"
        elif config.no_good_cuts:
            method = "no-good"
        elif config.intersection_cuts:
            method = "intersection"
        else:
            method = "branchandbound"

        average_cut_count = (
            np.mean(tracker.node_cut_count) if tracker.node_cut_count else 0
        )
        v_print(1, f"Violations: {tracker.violations}")
        v_print(1, f"Cut Count per Node: {tracker.node_cut_count}")
        v_print(1, f"Initial Solution given to callback: {config.solution}")
        v_print(1, f"Violations of initial solution: {tracker.solution_violations}")
        metadata = config.instance_metadata
        objective_multiplier = metadata.get("objective_multiplier", 1.0)
        raw_incumbent = objective
        raw_bound = best_bound
        recovered_incumbent = raw_incumbent * objective_multiplier
        recovered_bound = raw_bound * objective_multiplier
        result = {
            "problem": metadata.get("problem") or config.instance_type,
            "instance": os.path.splitext(os.path.basename(config.instance_file))[0],
            "n": metadata.get("n", ""),
            "density": metadata.get("density", ""),
            "replicate": metadata.get("replicate", ""),
            "formulation": metadata.get("formulation")
            or ("projected" if config.projected else "unprojected"),
            "linearcons": metadata.get("linearcons", ""),
            "method": method,
            "projected": int(config.projected),
            "status": status,
            "solver_status": solver_status,
            "has_incumbent": int(has_incumbent),
            "incumbent": recovered_incumbent,
            "bound": recovered_bound,
            "gap": gap,
            "gap_percent": gap * 100.0,
            "solve_time": tracker.solving_time,
            "bound_time": config.bound_time,
            "model_time_limit": config.time_lim,
            "total_time": tracker.solving_time + config.bound_time,
            "nodes_e": cpx.CPXgetnodecnt(g.env, g.m),
            "cut_num": tracker.cut_count,
            "cut_time": tracker.callback_time,
            "callback_count": tracker.callback_count,
            "cone_coefficient_time": tracker.cone_coefficient_time,
            "lower_level_constraints": self.n_lower_level_constraints,
            "objective_scale": metadata.get("objective_scale", 1.0),
            "raw_incumbent": raw_incumbent,
            "raw_bound": raw_bound,
        }
        print("result ,", ",".join(str(value) for value in result.values()))
        return result

    def solve(self):
        """Solve the knapsack problem."""
        # Set parameters
        config.set_cplex_parameters(g.env)

        try:
            self.set_callbacks()

            self.set_priorities()

            # Solve problem
            tracker.start_elapsed_time = cpx.CPXgettime(g.env)

            first_time = cpx.CPXgettime(g.env)
            try:
                cpx.CPXmipopt(g.env, g.m)
            finally:
                tracker.solving_time = cpx.CPXgettime(g.env) - first_time

            # Print solution
            return self.print_solution()

        except (CplexError, RuntimeError) as e:
            print(f"Error solving model: {e}")
            raise
        finally:
            # Always release the CPLEX resources, including error paths.
            cpx.CPXfreeprob(g.env, g.m)
            cpx.CPXcloseCPLEX(g.env)


#############################################
# Callback Classes
#############################################


class CutCallback(cpx.CutCallback):
    """
    Cut Callback class for CPLEX optimization.
    """

    def __init__(self, model):
        """Initialize the cut callback with the original model data."""
        super().__init__()
        # Initialize callback data from original model
        self.n = model.n
        self.p = model.p
        self.q = model.q
        self.r = model.r

        self.I = model.I
        self.J = model.J
        self.K = model.K

        self.b = model.b
        self.D = model.D
        self.C = model.C
        self.a = model.a
        self.beta = model.beta
        self.alpha = model.alpha
        self.ub = model.ub
        self.lb = model.lb
        self.ub_init = self.ub.copy()
        self.lb_init = self.lb.copy()
        
        self.bilinearities = model.bilinearities

    def __call__(self):
        """Callback function to be called by CPLEX during optimization."""
        tracker.start_callback_time = time.time()
        # Get time information
        self.get_time_information()

        try:
            tracker.callback_count += 1
            v_print(
                3,
                f"\n%%%%%%%%%%%%%%%%% Cut Callback Count {tracker.callback_count} %%%%%%%%%%%%%%%%%%%%%%%\n",
            )
            v_print(1, "Python callback called from event " + str(self.wherefrom))

            # Check time limit
            if self.check_time_limit():
                return cpx.CPX_CALLBACK_DEFAULT

            # Get node's problem
            self.lp = cpx.CPXgetcallbacknodelp(g.env, self.cbdata, self.wherefrom)

            # Write node problem to file for debugging
            if config.write_lps:
                cpx.CPXwriteprob(
                    g.env, self.lp, f"node_problem{tracker.callback_count}.lp"
                )

            # Check if node has already received too many cuts
            if not self.check_number_of_cuts():
                return cpx.CPX_CALLBACK_DEFAULT

            # Get node solution
            self.get_node_solution()

            global obj_val
            # Solve subproblem
            self.u_hat, obj_val = self.solve_subproblem()

            # Check feasibility of u_hat
            if self.check_feasibility():
                return cpx.CPX_CALLBACK_DEFAULT

            # If u_hat is not feasible, add a cut
            self.add_cut()

            # Update time information
            global callback_time
            tracker.callback_time += time.time() - tracker.start_callback_time
            v_print(
                2,
                f"Callback time this call: {time.time() - tracker.start_callback_time:.4f}s",
            )
            v_print(2, f"Total callback time so far: {tracker.callback_time:.4f}s")
            global total_time
            tracker.total_time += time.time() - self.start_time
            self.el_time = cpx.CPXgettime(g.env) - tracker.start_elapsed_time
            v_print(2, f"Total callback time: {tracker.total_time:.4f}s")
            v_print(2, f"Elapsed time from CPLEX: {self.el_time:.4f}s")
            v_print(
                3,
                "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n",
            )
        except:
            v_print(3, sys.exc_info()[0])
            raise

        return cpx.CPX_CALLBACK_SET

    def check_time_limit(self):
        """Return whether the main model has exhausted its time limit."""
        if self.el_time >= self.tilim:
            print(
                f"Elapsed time {self.el_time} exceeded time limit {self.tilim}. Stopping callback."
            )
            return True
        return False

    def check_number_of_cuts(self):
        """Check if the current node has already received the maximum number of cuts."""
        # Get node hash
        self.hash = self.get_node_hash()

        if self.hash not in tracker.node_hashes:
            tracker.node_hashes.append(self.hash)
            tracker.node_cut_count.append(0)
            tracker.node_initial_violation.append(0.0)
        if config.only_root_node:
            if self.hash != tracker.node_hashes[0]:
                v_print(
                    3,
                    "%%%%%%%%%%%%%%%% Not at root node. Skipping cuts. %%%%%%%%%%%%%%%%",
                )
                return False

        # Check node LP status
        self.stat = cpx.CPXgetstat(g.env, self.lp)
        if self.stat != 1:
            v_print(3, f"Warning: Node LP status is {self.stat}, expected 1 (optimal).")
            return False

        # Check if node has already received too many cuts
        if (
            tracker.node_cut_count[tracker.node_hashes.index(self.hash)]
            >= config.max_cuts
        ):
            v_print(
                3,
                f"%%%%%%%%%%%%%%%% Node has already received {tracker.node_cut_count[tracker.node_hashes.index(self.hash)]} cuts. Skipping further cuts.",
            )
            return False
        v_print(
            3,
            f"Node has received {tracker.node_cut_count[tracker.node_hashes.index(self.hash)]} cuts so far.",
        )
        return True

    def get_node_solution(self):
        """Get the solution of the current node's LP relaxation."""
        self.n_cols = cpx.CPXgetnumcols(g.env, self.lp)
        self.x = cpx.CPXgetx(g.env, self.lp, 0, self.n_cols - 1)
        v_print(1, f"Node solution (x,y,eta) = {self.x}")

    def get_time_information(self):
        """Get time information for the callback."""
        self.start_time = time.time()
        self.el_time = cpx.CPXgettime(g.env) - tracker.start_elapsed_time
        v_print(2, f"Elapsed time since optimization start: {self.el_time:.4f}s")
        self.tilim = cpx.CPXgetdblparam(g.env, cpx.CPX_PARAM_TILIM)

    def check_feasibility(self):
        """Check the feasibility of the solution u_hat obtained from the subproblem."""
        # Check feasibility of u_hat
        if config.instance_type == "knapsack":
            if (
                sum(self.a[i] * self.x[self.J[i]] for i in range(self.q+self.r))
                + obj_val
                > self.beta
            ):
                # if not bilevel feasible:
                v_print(1, f"{sum(self.a[i] * self.x[self.p + i] for i in range(self.q+self.r))} - {obj_val} > {self.beta - config.tolerance}")
                v_print(1, "u_hat is not upper-level feasible.")
                return False
            else:
                v_print(1, "u_hat is upper-level feasible.")
                v_print(
                    1,
                    f"{sum(self.a[i] * self.x[self.J[i]] for i in range(self.q))}, {obj_val}, {self.beta}",
                )
                return True
        elif config.instance_type == "bobilib":
            if (
                sum(self.a[i] * self.x[self.p + i] for i in range(self.q+self.r))
                - obj_val
                > self.beta
            ):
                # if not bilevel feasible:
                v_print(1, f"sum({self.a} * {self.x[self.p:self.p + self.q+self.r]}) + {obj_val} > {self.beta}")
                v_print(1, f"{sum(self.a[i] * self.x[self.p + i] for i in range(self.q+self.r))} + {obj_val} > {self.beta}")
                v_print(1, "u_hat is not upper-level feasible.")
                return False
            else:
                v_print(1, "u_hat is upper-level feasible.")
                v_print(1, f"sum({self.a} * {self.x[self.p:self.p + self.q+self.r]}) + {obj_val} <= {self.beta}")
                v_print(1, f"{sum(self.a[i] * self.x[self.p + i] for i in range(self.q+self.r))} + {obj_val} <= {self.beta}")
                return True

    def initialize_subproblem(self):
        """Initialize the subproblem model based on the current node solution and original model data."""
        # Create an empty model to start with
        self.sub_env = cpx.CPXopenCPLEX()
        self.sub_m = cpx.CPXcreateprob(self.sub_env, "Sub-Problem Model")

        # Add decision variables (u variables)
        var_names = [f"u_{i}" for i in range(self.q+self.r)]
        var_types = "I" * len(var_names)  # Integer variables
        # Calculate objective coefficients
        self.obj_coeffs = [0.0] * (len(var_names))
        for i in range(self.q+self.r):
            if config.instance_type == "knapsack":
                self.obj_coeffs[i] = self.b[i] * self.x[self.J[i]]
            elif config.instance_type == "bobilib":
                self.obj_coeffs[i] = self.a[i]
        # Add bounds
        lb = []
        ub = []
        for i in range(self.p, self.p + len(var_names)):
            lb.append(self.lb_init[i])
            ub.append(self.ub_init[i])
        # Add variables
        cpx.CPXnewcols(
            self.sub_env,
            self.sub_m,
            len(var_names),
            self.obj_coeffs,
            lb,
            ub,
            var_types,
            var_names,
        )
        # Set objective sense
        if config.instance_type == "knapsack":
            cpx.CPXchgobjsen(self.sub_env, self.sub_m, cpx.CPX_MAX)
        else:
            cpx.CPXchgobjsen(self.sub_env, self.sub_m, cpx.CPX_MIN)
        # Add constraints
        for j in range(len(self.D)):
            constraint_indices = []
            constraint_coeffs = []
            constraint_senses = []
            for i in range(len(self.D[j])):
                constraint_indices.append(i)
                constraint_coeffs.append(self.D[j][i])
            if config.instance_type == "knapsack":
                constraint_senses.append("L")
                if config.lower_level == "general":
                    rhs = self.alpha[j] + sum(
                        self.C[j][i] * self.x[self.J[i]] for i in range(len(self.J))
                    )
                else:
                    rhs = self.alpha[j] - sum(
                    self.C[j][i] * self.x[i] for i in range(len(self.I))
                )
            elif config.instance_type == "bobilib":
                constraint_senses.append("G")
                rhs = self.alpha[j] - sum(
                    self.C[j][i] * self.x[i] for i in range(len(self.I))
                )
            cpx.CPXaddrows(
                self.sub_env,
                self.sub_m,
                0,
                1,
                len(constraint_coeffs),
                [rhs],
                constraint_senses,
                [0],
                constraint_indices,
                constraint_coeffs,
                None,
                ["lower_capacity"],
            )
        try:
            config.sub_time_lim = max(
                0.0,
                cpx.CPXgetdblparam(g.env, cpx.CPX_PARAM_TILIM)
                - (cpx.CPXgettime(g.env) - tracker.start_elapsed_time),
            )
            config.set_cplex_parameters(self.sub_env)
        except (CplexError, RuntimeError) as e:
            v_print(3, f"Error setting subproblem time limit: {e}")

        if config.write_lps:
            cpx.CPXwriteprob(
                self.sub_env, self.sub_m, f"subproblem{tracker.callback_count}.lp"
            )
            v_print(
                1, f"Subproblem LP written to subproblem{tracker.callback_count}.lp"
            )

    def solve_subproblem(self):
        """Solve the subproblem and return the solution and objective value."""
        v_print(1, "\n********* Solving subproblem *********\n")

        # Initialize subproblem
        self.initialize_subproblem()

        subproblem_start = time.time()
        cpx.CPXmipopt(self.sub_env, self.sub_m)
        subproblem_time = time.time() - subproblem_start

        v_print(1, "Subproblem solved.")
        try:
            u_hat = cpx.CPXgetx(self.sub_env, self.sub_m, 0, self.q + self.r - 1)
            v_print(1, f"Subproblem solution (u_hat): {u_hat}")
            # Compute objective value manually from solution
            obj_val = sum(self.obj_coeffs[i] * u_hat[i] for i in range(self.q + self.r))
            if config.instance_type == "bobilib":
                obj_val = obj_val
            v_print(1, f"Subproblem objective coefficients: {self.obj_coeffs}")
            v_print(1, f"Subproblem objective value: {obj_val}\n")
        except Exception as e:
            v_print(3, f"Error retrieving subproblem solution: {e}")

        # Free model
        cpx.CPXfreeprob(self.sub_env, self.sub_m)

        # Free environment
        cpx.CPXcloseCPLEX(self.sub_env)

        v_print(2, f"Subproblem solve time: {subproblem_time:.4f}s")

        return u_hat, obj_val


class IntersectionCutCallback(CutCallback):
    """
    Intersection cut callback class for CPLEX optimization.
    """

    def __init__(self, model):
        """Initialize the intersection cut callback with the original model data."""
        super().__init__(model)
        if config.projected:
            self.eta_idx = model.eta_idx

    def add_cut(self):
        """Add an intersection cut to the current node's problem."""
        self.start_time = time.time()

        v_print(1, "\n********* Getting problem data *********\n")
        data_start = time.time()
        self.get_problem_data()
        data_time = time.time() - data_start
        v_print(2, f"Get problem data time: {data_time:.4f}s")

        v_print(1, "\n********* Getting rays *********\n")
        rays_start = time.time()
        self.rays = self.get_rays()
        rays_time = time.time() - rays_start
        v_print(2, f"Get rays time: {rays_time:.4f}s")

        v_print(1, "\n********* Getting robust free set *********\n")
        free_set_start = time.time()
        if self.bilinearities:
            self.F, self.g = self.rfs_with_bilinearities()
        else:
            self.F, self.g = self.rfs_without_bilinearities()
        free_set_time = time.time() - free_set_start
        v_print(2, f"Get robust free set time: {free_set_time:.4f}s")

        v_print(1, "\n********* Getting cone coefficients *********\n")
        cone_start = time.time()
        self.lam_num, self.lam_den = (
            self.get_cone_coefficients()
        )  # Now returns two arrays
        cone_time = time.time() - cone_start
        v_print(2, f"Get cone coefficients time: {cone_time:.4f}s")

        v_print(1, "\n********* Adding intersection cut *********\n")
        cut_start = time.time()
        self.add_intersection_cut()
        cut_time = time.time() - cut_start
        v_print(2, f"Add intersection cut time: {cut_time:.4f}s")

        total_time = time.time() - self.start_time
        v_print(2, f"Total intersection cut callback time: {total_time:.4f}s")

    def get_problem_data(self):
        """Retrieve problem data from the current node."""

        # Get number of columns and constraints
        self.n_rows = cpx.CPXgetnumrows(g.env, self.lp)
        v_print(1, f"self.n_cols = {self.n_cols}")
        v_print(1, f"self.n_rows = {self.n_rows}")

        # Get current solution
        self.x = cpx.CPXgetx(g.env, self.lp, 0, self.n_cols - 1)
        v_print(1, f"Node solution = {self.x}")

        # Get current bounds
        self.lb = cpx.CPXgetlb(g.env, self.lp, 0, self.n_cols - 1)
        self.ub = cpx.CPXgetub(g.env, self.lp, 0, self.n_cols - 1)
        v_print(1, f"LB = {self.lb}")
        v_print(1, f"UB = {self.ub}")

        # Get row sense
        self.sense = cpx.CPXgetsense(g.env, self.lp, 0, self.n_rows - 1)
        v_print(1, f"sense = {self.sense}")

        # Get basis
        (self.cstat, self.rstat) = cpx.CPXgetbase(g.env, self.lp)
        v_print(1, f"rstat = {self.rstat}")
        v_print(1, f"cstat = {self.cstat}")

        # Get basis head
        (self.head, self.x_B) = cpx.CPXgetbhead(g.env, self.lp)
        v_print(1, f"head = {self.head}")
        v_print(1, f"x_B = {self.x_B}")

    def get_rays(self):
        """Compute the rays of the current node's LP relaxation based on the basis information."""
        # Get non-basic variable indices
        self.non_basic_variable_indices = []

        # Non-basic structural variables (where cstat[j] != 1)
        self.non_basic_variable_indices.extend(
            [j for j in range(self.n_cols) if self.cstat[j] != 1]
        )

        # Non-basic slack variables (where rstat[i] != 1)
        self.non_basic_variable_indices.extend(
            [-i - 1 for i in range(self.n_rows) if self.rstat[i] != 1]
        )

        # Number of (non-)basic variables
        self.n_non_basic_variables = len(self.non_basic_variable_indices)

        # Compute the rays
        ray_start = time.time()
        v_print(1, f"NonBasic variables: {self.non_basic_variable_indices}")
        v_print(1, f"Basic variables: {self.head}\n")

        # Use NumPy for better performance
        rays = np.zeros((self.n_non_basic_variables, self.n_non_basic_variables))
        v_print(2, f"Size of rays matrix: {rays.shape}")

        # Pre-compute head indices for faster lookup
        head_dict = {val: idx for idx, val in enumerate(self.head)}

        for k in range(self.n_non_basic_variables):
            k_variable = self.non_basic_variable_indices[k]
            BinvAk = cpx.CPXbinvacol(g.env, self.lp, k_variable)
            v_print(1, f"BinvAk for variable {k_variable}: {BinvAk}")

            if k_variable >= 0:  # Structural variables
                rays[k, k_variable] = -1.0 if self.cstat[k_variable] == 2 else 1.0
                coeff = 1.0 if self.cstat[k_variable] == 2 else -1.0
            else:  # Slack Variables
                coeff = -1.0 if self.sense[-k_variable - 1] == "L" else 1.0

            for i in range(self.n_non_basic_variables):
                if i in head_dict:
                    rays[k][i] = coeff * BinvAk[head_dict[i]]
        
        #if config.instance_type == "bobilib":
        #    rays[:,self.eta_idx] = -rays[:,self.eta_idx]
        
        v_print(1, f"Final Rays = \n{np.array(rays)}")
        ray_time = time.time() - ray_start
        v_print(2, f"Ray computation time: {ray_time:.4f}s")
        return rays

    def rfs_with_bilinearities(self):
        """Compute the robust-free set with bilinearities."""
        # Fill matrix F
        F = [[0.0] * (len(self.x)) for _ in range(len(self.C) + 1)]
        # Fill first row
        F[0][-1] = -1.0

        for j in range(len(self.b)):
            F[0][self.J[j]] = self.b[j] * self.u_hat[j]

        for i in range(len(self.C)):
            for j in range(len(self.C[i])):
                F[i + 1][j] = -self.C[i][j]
        v_print(1, f"F = {F}")

        # Fill vector g
        g = [0.0] * (len(self.D) + 1)
        g[0] = 0.0

        # Pre-compute D * u_hat once
        D_u_hat = [
            sum(self.D[i][j] * self.u_hat[j] for j in self.I)
            for i in range(len(self.D))
        ]

        for i in range(len(self.D)):
            g[i + 1] = -self.alpha[i] - 1 + D_u_hat[i]

        v_print(1, f"g = {g}\n")

        return F, g

    def rfs_without_bilinearities(self):
        """Compute the robust-free set without bilinearities."""
        # Fill matrix F
        F = [[0.0] * (len(self.x)) for _ in range(len(self.C) + 1)]
        # Fill first row
        F[0][-1] = 1.0

        for i in range(len(self.C)):
            for j in self.I:
                F[i + 1][self.I[j]] = self.C[i][j]
        v_print(1, f"F = {F}")
    
        # Fill vector g
        g = [0.0] * (len(self.D) + 1)
        g[0] = sum(
            self.a[j] * self.u_hat[j]
            for j in range(self.r)
            )

        # Pre-compute D * u_hat once
        D_u_hat = [
            sum(self.D[i][j] * self.u_hat[j] for j in range(len(self.D[i])))
            for i in range(len(self.D))
        ]

        for i in range(len(self.D)):
            g[i + 1] = -(-self.alpha[i] + 1 + D_u_hat[i])

        v_print(1, f"g = {g}\n")

        return F, g
    
    def get_cone_coefficients(self):
        """Compute the cone coefficients lambda of the intersection cut based on the rays and the robust-free set."""
        cone_start = time.time()

        # Convert to NumPy arrays for efficient computation
        F_array = np.array(self.F)
        rays_array = np.array(self.rays)
        g_array = np.array(self.g)

        # Extract relevant part of x for F*x computation
        n_vars = F_array.shape[1]
        x_relevant = np.array(self.x[:n_vars])

        # Compute F*x using matrix-vector multiplication
        Fx_start = time.time()
        Fx = F_array @ x_relevant
        v_print(1, f"F_array = \n{F_array}\n")
        v_print(1, f"x_relevant = {x_relevant}\n")
        v_print(1, f"Fx = {Fx}\n")
        Fx_time = time.time() - Fx_start
        v_print(2, f"F*x computation time: {Fx_time:.4f}s")

        # Compute numerators: g - F*x
        numerators = g_array - Fx

        # Compute F*rays
        F_rays_start = time.time()

        # Create a mapping from variable indices to F columns
        F_relevant = np.zeros((F_array.shape[0], len(self.non_basic_variable_indices)))

        for k, var_idx in enumerate(self.non_basic_variable_indices):
            if var_idx >= 0 and var_idx < F_array.shape[1]:
                F_relevant[:, k] = F_array[:, var_idx]

        # Compute F_relevant * rays^T
        F_rays = F_array @ rays_array.T

        F_rays_time = time.time() - F_rays_start
        v_print(2, f"F*rays matrix computation time: {F_rays_time:.4f}s")
        v_print(1, f"F_array = \n{F_array}\n")
        v_print(1, f"F_relevant = \n{F_relevant}\n")
        v_print(1, f"Fx = {Fx}\n")
        v_print(1, f"F*rays = \n{F_rays}\n")

        # Store lambda as separate numerator and denominator arrays
        self.lam_num = []
        self.lam_den = []
        self.lam_indices = []  # Track which row index gives the minimum ratio

        lambda_start = time.time()
        epsilon = config.tolerance

        for k in range(len(self.rays)):
            # Find negative entries in F_rays[:, k]
            negative_mask = F_rays[:, k] < -epsilon

            if np.any(negative_mask):
                # Get the numerators and denominators for this ray
                nums = numerators[negative_mask]
                dens = F_rays[negative_mask, k]
                v_print(1, f"Ray {k}: Numerators = {nums}, Denominators = {dens}")

                # Find the minimum ratio (which gives lambda)
                ratios = nums / dens
                min_idx = np.argmin(ratios)

                # Get the original row index (before masking)
                original_indices = np.where(negative_mask)[0]
                chosen_row_idx = original_indices[min_idx]

                # Store the numerator and denominator separately
                self.lam_num.append(float(nums[min_idx]))
                self.lam_den.append(float(dens[min_idx]))
                self.lam_indices.append(int(chosen_row_idx))

                v_print(
                    1,
                    f"Ray {k}: row index = {chosen_row_idx}, lambda = {nums[min_idx]}/{dens[min_idx]} = {ratios[min_idx]:.6f}",
                )
            else:
                # No negative entries, lambda is infinite
                self.lam_num.append(float("inf"))
                self.lam_den.append(1.0)
                self.lam_indices.append(-1)  # -1 indicates no valid index

        lambda_time = time.time() - lambda_start
        v_print(2, f"Lambda computation time: {lambda_time:.4f}s")

        v_print(1, f"\nLambda numerators = {self.lam_num}")
        v_print(1, f"Lambda denominators = {self.lam_den}")
        v_print(1, f"Lambda row indices = {self.lam_indices}\n")

        cone_time = time.time() - cone_start
        v_print(2, f"Total cone coefficients computation time: {cone_time:.4f}s")

        global cone_coefficient_time
        tracker.cone_coefficient_time += cone_time
        v_print(
            2,
            f"Accumulated cone coefficient computation time: {tracker.cone_coefficient_time:.4f}s",
        )

        return self.lam_num, self.lam_den

    def compute_cut_coefficients(self):
        """Compute the coefficients of the intersection cut based on the rays, the cone coefficients lambda, and the problem data."""
        self.cut_indices = list(range(self.n_non_basic_variables))
        self.cut_coeffs = [0.0] * self.n_non_basic_variables
        self.rhs = 0.0

        rmatbeg, rmatind, rmatval = cpx.CPXgetrows(g.env, self.lp, 0, self.n_rows - 1)
        rmatrhs = cpx.CPXgetrhs(g.env, self.lp, 0, self.n_rows - 1)

        v_print(1, f"Rows: {rmatbeg, rmatind, rmatval}")
        v_print(1, f"RHS: {rmatrhs}")

        min_lambda_num = config.tolerance  # Minimum lambda numerator to avoid issues

        for i in range(len(self.non_basic_variable_indices)):
            if abs(self.lam_num[i]) < min_lambda_num:
                v_print(
                    1,
                    f"Warning: Lambda[{i}] = {self.lam_num[i]}/{self.lam_den[i]} is invalid, skipping",
                )
                continue

            # Structural variables
            if self.non_basic_variable_indices[i] >= 0:
                v_print(
                    1,
                    f"\nStructural variable: {self.non_basic_variable_indices[i]}, Lambda: {self.lam_num[i]}/{self.lam_den[i]}, Bound: {self.cstat[self.non_basic_variable_indices[i]]}",
                )

                if (
                    self.cstat[self.non_basic_variable_indices[i]] == 2
                ):  # At upper bound
                    # Add: (1/lambda) * coeff = (lam_den[i] / lam_num[i]) * coeff
                    self.cut_coeffs[self.non_basic_variable_indices[i]] += (
                        self.lam_den[i] / self.lam_num[i]
                    )
                    # RHS += (1/lambda) * bound = (lam_den[i] / lam_num[i]) * ub
                    self.rhs += (
                        self.lam_den[i] * self.ub[self.non_basic_variable_indices[i]]
                    ) / self.lam_num[i]
                else:  # At lower bound
                    # Subtract: (1/lambda) * coeff
                    self.cut_coeffs[self.non_basic_variable_indices[i]] -= (
                        self.lam_den[i] / self.lam_num[i]
                    )
                    # RHS -= (1/lambda) * bound
                    self.rhs -= (
                        self.lam_den[i] * self.lb[self.non_basic_variable_indices[i]]
                    ) / self.lam_num[i]

            # Slack variables
            elif self.non_basic_variable_indices[i] < 0:
                v_print(
                    1,
                    f"\nSlack variable: {self.non_basic_variable_indices[i]}, Lambda: {self.lam_num[i]}/{self.lam_den[i]}, Sense: {self.sense[-self.non_basic_variable_indices[i]-1]}",
                )
                start_pos = rmatbeg[-self.non_basic_variable_indices[i] - 1]
                if -self.non_basic_variable_indices[i] - 1 < len(rmatbeg) - 1:
                    end_pos = rmatbeg[-self.non_basic_variable_indices[i]]
                else:
                    end_pos = len(rmatval)
                row_length = end_pos - start_pos

                if self.sense[-self.non_basic_variable_indices[i] - 1] == "L":
                    for j in range(row_length):
                        # coeff += (1/lambda) * rmatval = (lam_den[i] / lam_num[i]) * rmatval
                        self.cut_coeffs[rmatind[j + start_pos]] += (
                            self.lam_den[i] * rmatval[j + start_pos]
                        ) / self.lam_num[i]
                    # RHS += (1/lambda) * rhs = (lam_den[i] / lam_num[i]) * rhs
                    self.rhs += (
                        self.lam_den[i]
                        * rmatrhs[-self.non_basic_variable_indices[i] - 1]
                    ) / self.lam_num[i]
                else:
                    for j in range(row_length):
                        # coeff -= (1/lambda) * rmatval
                        self.cut_coeffs[rmatind[j + start_pos]] -= (
                            self.lam_den[i] * rmatval[j + start_pos]
                        ) / self.lam_num[i]
                    # RHS -= (1/lambda) * rhs
                    self.rhs -= (
                        self.lam_den[i]
                        * rmatrhs[-self.non_basic_variable_indices[i] - 1]
                    ) / self.lam_num[i]

        self.rhs += -1

        self.cut_coeffs = [float(coef) for coef in self.cut_coeffs]

        v_print(1, f"\nFinal Cut Indices: {self.cut_indices}")
        v_print(1, f"Final Cut Coefficients: {self.cut_coeffs}")
        v_print(1, f"Final RHS: {self.rhs}\n")

    def check_violation(self):
        """Check the violation of the intersection cut and decide whether to add it."""
        v_print(
            1,
            f"Checking violation: |{self.cut_coeffs} * {self.x} - {self.rhs}| / ||{self.cut_coeffs}||_2 > {config.tolerance}?",
        )
        v_print(
            1,
            f"Checking violation: {abs(sum(self.cut_coeffs[i] * self.x[i] for i in range(len(self.cut_coeffs)))- self.rhs)} / {np.linalg.norm(self.cut_coeffs)} > {config.tolerance}?",
        )
        non_zero_coeffs = [
            abs(coeff) for coeff in self.cut_coeffs if abs(coeff) > config.tolerance
        ]
        # Compute the ratio of the maximum to minimum absolute non-zero coefficient
        if non_zero_coeffs:
            ratio = max(non_zero_coeffs) / min(non_zero_coeffs)
            v_print(1, f"Coefficient ratio (max/min absolute values): {ratio}")
            current_violation = float(
                abs(
                    sum(self.cut_coeffs[i] * self.x[i] for i in range(len(self.cut_coeffs)))
                    - self.rhs
                )
                / np.linalg.norm(self.cut_coeffs)
            )
            if tracker.node_initial_violation[tracker.node_hashes.index(self.hash)] == 0.0:
                tracker.node_initial_violation[tracker.node_hashes.index(self.hash)] = (
                    current_violation
                )
                v_print(
                    1, f"Setting initial violation for this node to {current_violation}"
                )
            v_print(1, f"Current violation: {current_violation}")
            v_print(1, f"Node solution: {self.x}")
            tracker.violations.append(current_violation)
        else:
            current_violation = 0.0
            ratio = 1.0  # No non-zero coefficients
            v_print(1, "All coefficients are zero")
            
        
        # Compute the violation of the initial solution if provided
        if config.solution:
            v_print(1, f"Initial Solution given to callback: {config.solution}")
            initial_solution_violation = float(
                sum(
                    self.cut_coeffs[i] * config.solution[i]
                    for i in range(len(self.cut_coeffs))
                )
                - self.rhs
            )
            tracker.solution_violations.append(initial_solution_violation)
            v_print(1, f"Initial solution violation: {initial_solution_violation}")

        min_violation = config.tolerance * (abs(self.rhs) + 1)
        relative_violation = 0.1

        # Only add the cut if the violation is significant relative to the initial violation 
        # and above an absolute threshold, and if the coefficient ratio is not too large
        return (
            current_violation
            >= relative_violation
            * tracker.node_initial_violation[tracker.node_hashes.index(self.hash)]
            and current_violation >= min_violation
            and ratio < 1e6
        )

    def add_intersection_cut(self):
        """Add intersection cut"""
        self.compute_cut_coefficients()

        tmp_sense = "L"

        # Check if the cut should be added
        if not self.check_violation():
            v_print(3, f"Cut not added as it is not violated enough.")
            return cpx.CPX_CALLBACK_DEFAULT

        # Add the cut using CPLEX
        cpx.CPXcutcallbackaddlocal(
            g.env,
            self.cbdata,
            self.wherefrom,
            len(self.cut_indices),
            self.rhs,
            tmp_sense,
            self.cut_indices,
            self.cut_coeffs,
        )

        # Update cut counts
        tracker.node_cut_count[tracker.node_hashes.index(self.hash)] += 1
        tracker.cut_count += 1
        v_print(3, f"Cut count updated: {tracker.cut_count}")
        return cpx.CPX_CALLBACK_DEFAULT


class InterdictionCutCallback(CutCallback):
    """
    Interdiction cut callback class for CPLEX optimization.
    """

    def __init__(self, model):
        """Initialize the interdiction cut callback with the original model data."""
        super().__init__(model)

    def add_cut(self):
        """Generate and add an interdiction cut to the main model."""
        v_print(1, "Generating interdiction cut...")

        # Calculate the coefficients for the interdiction cut
        cut_indices = []
        cut_coeffs = []

        # Add coefficients for y variables: bi * ui for each i
        for i in range(self.q):
            cut_indices.append(self.J[i])  # y variable index
            coeff = self.b[i] * self.u_hat[i]
            cut_coeffs.append(coeff)

        # Add coefficients for z variables (McCormick): -bk * uk for each k
        for k in range(min(self.p, self.q)):
            z_idx = self.n + 1 + k  # z_k variable index
            cut_indices.append(z_idx)
            coeff = -self.b[k] * self.u_hat[k]
            cut_coeffs.append(coeff)

        # Add coefficient for η (auxiliary variable): -1
        cut_indices.append(self.n)  # eta variable index (stored as self.n)
        cut_coeffs.append(-1.0)

        # Right-hand side is 0 (since we rearranged to ≤ 0 form)
        rhs = 0.0
        sense = "L"  # Less than or equal constraint
        nzcnt = len(cut_indices)  # Number of non-zero coefficients

        v_print(1, f"Interdiction cut constraint:")
        v_print(1, f"  Variables: {len(cut_indices)} terms")
        v_print(1, f"  Y coefficients: {cut_coeffs[:self.q]}")
        if len(cut_coeffs) > self.q:
            v_print(1, f"  Z coefficients: {cut_coeffs[self.q:-1]}")
        v_print(1, f"  Eta coefficient: {cut_coeffs[-1]}")
        v_print(1, f"  RHS: <= {rhs}")

        # Add the constraint using CPLEX to the main model
        try:
            # Add the lazy constraint to the main model
            v_print(1, "Adding interdiction cut to main model...")

            cpx.CPXcutcallbackadd(
                g.env,
                self.cbdata,
                self.wherefrom,
                nzcnt,
                rhs,
                sense,
                cut_indices,
                cut_coeffs,
                0,
            )

            # Write the node model to an LP file for inspection
            self.lp = cpx.CPXgetcallbacknodelp(g.env, self.cbdata, self.wherefrom)

            # Increment the interdiction cut count
            tracker.cut_count += 1

            v_print(
                1,
                f"Successfully added interdiction cut #{tracker.cut_count} to main model",
            )

        except Exception as e:
            v_print(3, f"Error adding interdiction cut: {e}")

        v_print(1, "Interdiction cut callback completed")


class NoGoodCutCallback(CutCallback):
    """
    No-good cut callback class for CPLEX optimization.
    """

    def __init__(self, model):
        """Initialize the no-good cut callback with the original model data."""
        super().__init__(model)

    def add_cut(self):
        """Generate and add a no-good cut to the main model."""

        v_print(1, "Generating no-good cut...")

        # Calculate the coefficients for the no-good cut
        cut_indices = []
        cut_coeffs = []
        rhs = 1.0
        sense = "G"

        for i in range(self.n):
            cut_indices.append(i)
            if abs(self.x[i] - 1.0) < config.tolerance:
                cut_coeffs.append(-1.0)
                rhs -= 1.0
            else:
                cut_coeffs.append(1.0)

        v_print(1, f"No-good cut constraint:")
        v_print(1, f"  Variables: {len(cut_indices)} terms")
        v_print(1, f"  Coefficients: {cut_coeffs}")
        v_print(1, f"  RHS: <= {rhs}")

        # Add the constraint using CPLEX to the main model
        try:
            # Add the lazy constraint to the main model
            v_print(1, "Adding no-good cut to main model...")

            cpx.CPXcutcallbackaddlocal(
                g.env,
                self.cbdata,
                self.wherefrom,
                len(cut_indices),
                rhs,
                sense,
                cut_indices,
                cut_coeffs,
            )

            # Update cut counts
            tracker.node_cut_count[tracker.node_hashes.index(self.hash)] += 1
            tracker.cut_count += 1
            v_print(3, f"Cut count updated: {tracker.cut_count}")
            v_print(1, "Successfully added no-good cut to main model")

        except Exception as e:
            v_print(3, f"Error adding no-good cut: {e}")

        v_print(1, "No-good cut callback completed")


class MyIncumbentCallback(cpx.IncumbentCallback):
    """
    Incumbent callback class for CPLEX optimization.
    """

    def __init__(self, model):
        """Initialize the incumbent callback with the original model data."""
        super().__init__()

        # Initialize callback data from original model
        self.n = model.n
        self.p = model.p
        self.q = model.q
        self.r = model.r

        self.I = model.I
        self.J = model.J
        self.K = model.K

        self.b = model.b
        self.D = model.D
        self.C = model.C
        self.a = model.a
        self.beta = model.beta
        self.alpha = model.alpha
        self.ub = model.ub
        self.lb = model.lb
        self.ub_init = self.ub.copy()
        self.lb_init = self.lb.copy()
        
        self.bilinearities = model.bilinearities

    def __call__(self):
        """Callback function that is called by CPLEX when a new incumbent solution is found."""
        tracker.callback_count += 1
        callback_start = time.time()
        v_print(3, "*************************************")
        v_print(3, "      IncumbentCallback Begin        ")
        v_print(3, "*************************************")

        v_print(1, "Python callback called from event " + str(self.wherefrom))
        self.el_time = cpx.CPXgettime(g.env) - tracker.start_elapsed_time
        v_print(2, f"Elapsed time since optimization start: {self.el_time:.4f}s")
        self.tilim = cpx.CPXgetdblparam(g.env, cpx.CPX_PARAM_TILIM)

        if self.check_time_limit():
            self.isfeas = False
            return cpx.CPX_CALLBACK_SET

        # Get Node's Problem
        self.lp = cpx.CPXgetcallbacknodelp(g.env, self.cbdata, self.wherefrom)
        if config.write_lps:
            lp_filename = f"nodeproblem_{tracker.callback_count}.lp"
            cpx.CPXwriteprob(g.env, self.lp, lp_filename, None)
            v_print(1, f"Node LP written to {lp_filename} for inspection.")

        self.n_cols = cpx.CPXgetnumcols(g.env, self.lp)
        self.x = self.get_x(self.lp)
        v_print(1, f"Upper-level solution x: {self.x}")

        global obj_val
        v_print(1, f"Subproblem printed to LP file {tracker.callback_count} for inspection.")
        self.u_hat, obj_val, sub_status = self.solve_subproblem()

        if sub_status != cpx.CPXMIP_OPTIMAL:
            v_print(
                1,
                f"Subproblem status {sub_status} is not optimal; rejecting incumbent candidate.",
            )
            self.isfeas = False
        elif not self.check_feasibility():
            self.isfeas = False
        else:
            self.isfeas = True

        callback_time = time.time() - callback_start
        v_print(2, f"Incumbent callback time: {callback_time:.4f}s")

        v_print(3, "*************************************")
        v_print(3, "      IncumbentCallback End          ")
        v_print(3, "*************************************")

        return cpx.CPX_CALLBACK_SET

    def check_time_limit(self):
        """Return whether the main model has exhausted its time limit."""
        if self.el_time >= self.tilim:
            print(
                f"Elapsed time {self.el_time} exceeded time limit {self.tilim}. Stopping callback."
            )
            return True
        return False

    def check_feasibility(self):
        """Check the feasibility of the current solution with respect to the upper-level constraints."""
        # Check feasibility of u_hat
        if config.instance_type == "knapsack":
            if (
                sum(self.a[i] * self.x[self.J[i]] for i in range(self.q+self.r))
                + obj_val
                > self.beta
            ):
                # if not bilevel feasible:
                v_print(1, f"{sum(self.a[i] * self.x[self.p + i] for i in range(self.q+self.r))} - {obj_val} > {self.beta - config.tolerance}")
                v_print(1, "u_hat is not upper-level feasible.")
                return False
            else:
                v_print(1, "u_hat is upper-level feasible.")
                v_print(
                    1,
                    f"{sum(self.a[i] * self.x[self.J[i]] for i in range(self.q))}, {obj_val}, {self.beta}",
                )
                return True
        elif config.instance_type == "bobilib":
            if (
                sum(self.a[i] * self.x[self.p + i] for i in range(self.q+self.r))
                - obj_val
                > self.beta
            ):
                # if not bilevel feasible:
                v_print(1, f"sum({self.a} * {self.x[self.p:self.p + self.q+self.r]}) + {obj_val} > {self.beta}")
                v_print(1, f"{sum(self.a[i] * self.x[self.p + i] for i in range(self.q+self.r))} + {obj_val} > {self.beta}")
                v_print(1, "u_hat is not upper-level feasible.")
                return False
            else:
                v_print(1, "u_hat is upper-level feasible.")
                v_print(1, f"sum({self.a} * {self.x[self.p:self.p + self.q+self.r]}) + {obj_val} <= {self.beta}")
                v_print(1, f"{sum(self.a[i] * self.x[self.p + i] for i in range(self.q+self.r))} + {obj_val} <= {self.beta}")
                return True

    def initialize_subproblem(self):
        """Initialize the subproblem model based on the current solution x and the original problem data."""
        # Create an empty model to start with
        self.sub_env = cpx.CPXopenCPLEX()
        self.sub_m = cpx.CPXcreateprob(self.sub_env, "Sub-Problem Model")

        # Add decision variables (u variables)
        var_names = [f"u_{i}" for i in range(self.q+self.r)]
        var_types = "I" * len(var_names)  # Integer variables
        # Calculate objective coefficients
        self.obj_coeffs = [0.0] * (len(var_names))
        for i in range(self.q+self.r):
            if config.instance_type == "knapsack":
                self.obj_coeffs[i] = self.b[i] * self.x[self.J[i]]
            elif config.instance_type == "bobilib":
                self.obj_coeffs[i] = self.a[i]
        # Add bounds
        lb = []
        ub = []
        for i in range(self.p, self.p + len(var_names)):
            lb.append(self.lb_init[i])
            ub.append(self.ub_init[i])
        # Add variables
        cpx.CPXnewcols(
            self.sub_env,
            self.sub_m,
            len(var_names),
            self.obj_coeffs,
            lb,
            ub,
            var_types,
            var_names,
        )
        # Set objective sense
        if config.instance_type == "knapsack":
            cpx.CPXchgobjsen(self.sub_env, self.sub_m, cpx.CPX_MAX)
        else:
            cpx.CPXchgobjsen(self.sub_env, self.sub_m, cpx.CPX_MIN)
        # Add constraints
        for j in range(len(self.D)):
            constraint_indices = []
            constraint_coeffs = []
            constraint_senses = []
            for i in range(len(self.D[j])):
                constraint_indices.append(i)
                constraint_coeffs.append(self.D[j][i])
            if config.instance_type == "knapsack":
                constraint_senses.append("L")
                if config.lower_level == "general":
                    rhs = self.alpha[j] + sum(
                        self.C[j][i] * self.x[self.J[i]] for i in range(len(self.J))
                    )
                else:
                    rhs = self.alpha[j] - sum(
                    self.C[j][i] * self.x[i] for i in range(len(self.I))
                )
            elif config.instance_type == "bobilib":
                constraint_senses.append("G")
                rhs = self.alpha[j] - sum(
                    self.C[j][i] * self.x[i] for i in range(len(self.I))
                )
            cpx.CPXaddrows(
                self.sub_env,
                self.sub_m,
                0,
                1,
                len(constraint_coeffs),
                [rhs],
                constraint_senses,
                [0],
                constraint_indices,
                constraint_coeffs,
                None,
                ["lower_capacity"],
            )

        config.sub_time_lim = max(
            0.0,
            cpx.CPXgetdblparam(g.env, cpx.CPX_PARAM_TILIM)
            - (cpx.CPXgettime(g.env) - tracker.start_elapsed_time),
        )
        config.set_cplex_parameters(self.sub_env)

        if config.write_lps:
            cpx.CPXwriteprob(
                self.sub_env, self.sub_m, f"subproblem{tracker.callback_count}.lp"
            )
            v_print(
                1, f"Subproblem LP written to subproblem{tracker.callback_count}.lp"
            )

    def solve_subproblem(self):
        """Solve the subproblem and return its solution, objective, and status."""
        v_print(1, "\n********* Solving subproblem *********\n")

        self.initialize_subproblem()

        subproblem_start = time.time()
        cpx.CPXmipopt(self.sub_env, self.sub_m)
        subproblem_time = time.time() - subproblem_start
        sub_status = cpx.CPXgetstat(self.sub_env, self.sub_m)

        v_print(1, "Subproblem solved.")
        u_hat = None
        obj_val = None
        if sub_status == cpx.CPXMIP_OPTIMAL:
            try:
                u_hat = cpx.CPXgetx(
                    self.sub_env, self.sub_m, 0, self.q + self.r - 1
                )
                v_print(1, f"Subproblem solution (u_hat): {u_hat}")
                # Compute objective value manually from solution
                obj_val = sum(
                    self.obj_coeffs[i] * u_hat[i]
                    for i in range(self.q + self.r)
                )
                v_print(1, f"Subproblem objective coefficients: {self.obj_coeffs}")
                v_print(1, f"Subproblem objective value: {obj_val}\n")
            except Exception as e:
                sub_status = None
                v_print(3, f"Error retrieving subproblem solution: {e}")

        # Free model
        cpx.CPXfreeprob(self.sub_env, self.sub_m)

        # Free environment
        cpx.CPXcloseCPLEX(self.sub_env)

        v_print(2, f"Subproblem solve time: {subproblem_time:.4f}s")

        return u_hat, obj_val, sub_status
