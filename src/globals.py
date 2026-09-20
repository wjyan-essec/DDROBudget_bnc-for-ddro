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


from python_to_cplex_c_api import cplex_c_api_wrapper as cpx

env = None
m = None

class Config:
    """Configuration class to hold parameters and settings for the optimization process."""
    def __init__(
        self,
        write_lps,
        tolerance,
        max_cuts,
        only_root_node,
        cplex_cuts,
        verbose_level,
        time_lim,
        bound_time=0.0,
        total_time_lim=None,
        instance_metadata=None,
        interdiction_cuts=True,
        intersection_cuts=True,
        nogood_cuts=True,
        branchandbound=False,
        projected=False,
        instance_type="general",
        instance_file="",
        lower_level="",
        separation="integer",
        solution=None,
    ):
        """Initialize the configuration with given parameters."""

        # Configuration parameters
        self.write_lps = write_lps
        self.tolerance = tolerance
        self.max_cuts = max_cuts
        self.only_root_node = only_root_node
        self.cplex_cuts = cplex_cuts
        self.verbose = verbose_level
        self.time_lim = time_lim
        self.bound_time = bound_time
        self.total_time_lim = time_lim if total_time_lim is None else total_time_lim
        self.instance_metadata = {} if instance_metadata is None else instance_metadata
        self.sub_time_lim = time_lim
        self.separation = separation

        # Cut parameters
        self.interdiction_cuts = interdiction_cuts
        self.intersection_cuts = intersection_cuts
        self.no_good_cuts = nogood_cuts
        self.branchandbound = branchandbound
        self.projected = projected

        # Instance parameter
        self.instance_type = instance_type
        self.instance_file = instance_file
        self.lower_level = lower_level

        # Create verbose print function
        self._v_print = self._create_v_print()

        # Solution parameter
        self.solution = solution

    def _create_v_print(self):
        """Create a verbose print function based on the verbosity level."""
        def _v_print(*verb_args):
            if verb_args[0] > (3 - self.verbose):
                print(verb_args[1])

        return _v_print

    def set_cplex_parameters(self, env):
        """Set CPLEX parameters based on the configuration."""
        # Set CPLEX parameters based on configuration
        int_params = [
            (cpx.CPX_PARAM_SCRIND, cpx.CPX_ON),  # Screen indicator
            (cpx.CPX_PARAM_MIPDISPLAY, 5),  # MIP display level
            (cpx.CPX_PARAM_MIPCBREDLP, cpx.CPX_OFF),  # Reduce LP in MIP callbacks
            (cpx.CPX_PARAM_PREIND, cpx.CPX_OFF),  # Disable presolve
            (cpx.CPX_PARAM_STARTALG, cpx.CPX_ALG_PRIMAL),  # Start algorithm
            (cpx.CPX_PARAM_SUBALG, cpx.CPX_ALG_DUAL),  # Sub-algorithm
            (cpx.CPX_PARAM_HEURFREQ, -1),  # Disable heuristics
            (cpx.CPX_PARAM_THREADS, 1),  # Single thread
            (cpx.CPX_PARAM_CUTPASS, self.cplex_cuts),  # No cut passes
            (cpx.CPX_PARAM_PRELINEAR, 0),  # Disable presolve linearization
            (cpx.CPX_PARAM_MIPORDIND, cpx.CPX_ON),  # MIP order indicator
            (cpx.CPX_PARAM_PROBE, 0),  # Probe
            (cpx.CPX_PARAM_REDUCE, 0),  # Reduce
            (cpx.CPX_PARAM_VARSEL, cpx.CPX_VARSEL_DEFAULT),  # Variable selection
            (cpx.CPX_PARAM_NUMERICALEMPHASIS, 1),  # Numerical emphasis
            (cpx.CPX_PARAM_CLOCKTYPE, 2),  # Use elapsed time
        ]
        if self.verbose == 0:
            int_params[0] = (cpx.CPX_PARAM_SCRIND, cpx.CPX_OFF)
            int_params[1] = (cpx.CPX_PARAM_MIPDISPLAY, 0)

        for param, value in int_params:
            cpx.CPXsetintparam(env, param, value)

        # Set double parameters
        cpx.CPXsetdblparam(
            env, cpx.CPX_PARAM_EPGAP, self.tolerance
        )  # MIP optimality gap
        cpx.CPXsetdblparam(env, cpx.CPX_PARAM_TILIM, self.sub_time_lim)  # Time limit
        cpx.CPXsetdblparam(
            env, cpx.CPX_PARAM_EPINT, self.tolerance
        )  # Integer feasibility tolerance


class Tracker:
    """Tracker class to keep track of various metrics during the optimization process."""
    def __init__(self):
        """Initialize the tracker with default values."""
        # Time tracking
        self.total_time = 0.0
        self.callback_time = 0.0
        self.cone_coefficient_time = 0.0
        self.start_elapsed_time = 0.0
        self.solving_time = 0.0
        self.start_callback_time = 0.0

        # Cut tracking
        self.cut_count = 0
        self.node_cut_count = []
        self.node_hashes = []
        self.node_initial_violation = []
        self.violations = []
        self.solution_violations = []

        self.interdiction_cut_count = 0
        self.intersection_cut_count = 0
        self.nogood_cut_count = 0

        # Callback tracking
        self.callback_count = 0

        # Subproblem tracking
        self.sub_count = 0
