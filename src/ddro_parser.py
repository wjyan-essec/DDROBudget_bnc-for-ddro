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

from globals import Config, Tracker
import globals as g

from mps_aux_parser import MPS_AUX_Parser
from parse_knapsack import *
import numpy as np


class Parser:
    """Base class for parsers."""

    def __init__(self, instance_file, config: Config):
        """Initialize the parser with the instance file and configuration."""
        self.instance_file = instance_file
        self.config = config
        
        # Initialize verbose print function
        global v_print
        v_print = config._v_print

class MPSAUXParser(Parser):
    """Parser for MPS and AUX files."""

    def __init__(self, instance_file, config: Config):
        """Initialize the MPSAUXParser with the instance file and configuration."""
        super().__init__(instance_file, config)

        self.aux_file = instance_file.replace(".mps", ".aux")

    def parse_bobilib_instance(self):
        """Parse the MPS and AUX files."""
        if self.config.instance_file.endswith(".mps"):
            aux_file = self.config.instance_file.replace(".mps", ".aux")
        elif self.config.instance_file.endswith(".mps.gz"):
            aux_file = self.config.instance_file.replace(".mps.gz", ".aux")
        mps_aux_reader = MPS_AUX_Parser(self.config.instance_file, aux_file)
        mps_aux_reader.read()
        instance_data = mps_aux_reader.get_instance_data()

        # Set dimensions
        self.p = instance_data.nr_ul_vars
        self.q = 0
        self.r = instance_data.nr_ll_vars
        self.n = self.p + self.r

        # Set variable and constraint names
        self.var_names = instance_data.ul_vars + instance_data.ll_vars
        self.ul_constr_names = instance_data.ul_constraints
        self.ll_constr_names = instance_data.ll_constraints

        # Set objective coefficients
        self.leader_obj_coeffs = list(instance_data.arrays["leader_ul_obj"])
        self.follower_obj_coeffs = list(instance_data.arrays["follower_ul_obj"])
        self.obj_coeffs = []
        for i in range(self.p):
            self.obj_coeffs.append(float(self.leader_obj_coeffs[i]))
        for i in range(self.q):
            self.obj_coeffs.append(float(self.follower_obj_coeffs[i]))

        # Set constraint matrices and right-hand sides
        self.leader_ul_constrs = instance_data.arrays["leader_ul_constr_mat"]
        self.follower_ul_constrs = instance_data.arrays["follower_ul_constr_mat"]
        self.leader_ll_constrs = instance_data.arrays["leader_ll_constr_mat"]
        self.follower_ll_constrs = instance_data.arrays["follower_ll_constr_mat"]
        
        self.ul_rhs = instance_data.arrays["ul_rhs"]
        self.ll_rhs = instance_data.arrays["ll_rhs"]

        # Set variable bounds
        self.leader_lbs = instance_data.arrays["leader_lbs"]
        self.leader_ubs = instance_data.arrays["leader_ubs"]
        self.follower_lbs = instance_data.arrays["follower_lbs"]
        self.follower_ubs = instance_data.arrays["follower_ubs"]
        
        self.ub = np.concatenate((self.leader_ubs, self.follower_ubs))
        self.lb = np.concatenate((self.leader_lbs, self.follower_lbs))
        
        # Set initial problem data
        self.c = []
        self.h = []
        for i in range(self.p):
            self.h.append(float(instance_data.arrays["leader_ul_obj"][i]))
        for i in range(self.r):
            self.c.append(float(instance_data.arrays["follower_ul_obj"][i]))

        self.beta = 0
        self.a = instance_data.arrays["follower_ll_obj"]
        self.b = instance_data.arrays["leader_ul_constr_mat"]
        self.n_lower_level_constraints = instance_data.nr_ll_constrs
        self.D = instance_data.arrays["follower_ll_constr_mat"]
        self.C = instance_data.arrays["leader_ll_constr_mat"]
        self.alpha = instance_data.arrays["ll_rhs"]
        
        # Fill A and gamma
        self.A_upper = np.hstack((self.leader_ul_constrs, self.follower_ul_constrs))
        self.A_lower = np.hstack((self.leader_ll_constrs, self.follower_ll_constrs))
        self.A = np.vstack((self.A_upper, self.A_lower))
        self.gamma = np.hstack((self.ul_rhs, self.ll_rhs))
        
        # Fill I and J
        self.I = range(self.p)
        self.J = []
        self.K = range(self.p, self.p + self.r)
        
        # Collect all problem data into a dictionary
        problem_data = {
            "c": self.c,
            "h": self.h,
            "beta": self.beta,
            "A": self.A,
            "gamma": self.gamma,
            "a": self.a,
            "b": self.b,
            "D": self.D,
            "C": self.C,
            "alpha": self.alpha,
            "I": self.I,
            "J": self.J,
            "K": self.K,
            "n": self.n,
            "p": self.p,
            "q": self.q,
            "r": self.r,
            "lb": self.lb,
            "ub": self.ub,
            "n_lower_level_constraints": self.n_lower_level_constraints,
            "instance_file": self.instance_file,
            "bilinearities": False,
            # Preserve the source names so application-specific adapters can
            # map semantic variables (for example, TSP edges) to master
            # columns without relying on a fragile positional convention.
            "var_names": self.var_names,
            "ul_var_names": list(instance_data.ul_vars),
            "ll_var_names": list(instance_data.ll_vars),
        }
        return problem_data

class KNAPSACKParser(Parser):
    """Parser for knapsack problem instances."""

    def __init__(self, instance_file, config: Config):
        """Initialize the KNAPSACKParser with the instance file and configuration."""
        super().__init__(instance_file, config)

    def parse_knapsack_instance(self):
        """Parse the knapsack instance from the file."""
        # Parse .kp file
        c, h, beta, a, b, D, alpha, C, instance_type, instance_file = parse_knapsack(
            self.instance_file
        )

        # Set initial problem data
        self.c = [-c_i for c_i in c]  # Negate c for minimization
        self.h = [float(h_i) for h_i in h]
        self.beta = beta
        self.a = a
        self.b = b
        self.n_lower_level_constraints = len(D)
        self.D = D
        self.C = C
        self.alpha = alpha

        # Add additional data for interdiction instances
        if self.config.lower_level == "interdiction":
            self.add_interdiction_data()

        v_print(1, f"D: {self.D}")
        v_print(1, f"C: {self.C}")
        v_print(1, f"alpha: {self.alpha}")

        # Set dimensions based on the number of leader and follower variables
        if self.config.lower_level == "interdiction":
            self.p = len(c)
            self.q = len(h)
            self.n = self.p + self.q  # Total number of variables (x_i and y_i)
        else:
            self.p = 0
            self.q = len(c)
            self.n = self.q  # Total number of variables (only y_i)
        
        # Set variable bounds (0-1 for knapsack variables)
        self.ub = [1.0] * self.n  # Upper bounds for x_i and y_i
        self.lb = [0.0] * self.n  # Lower bounds for x_i
        
        # Fill I and J based on the number of leader and follower variables
        self.I = list(range(self.p))
        self.J = [i + self.p for i in range(self.q)]
        self.K = []  # No non-affected variables in knapsack problem
        
        # Collect all problem data into a dictionary
        problem_data = {
            "c": self.c,
            "h": self.h,
            "beta": self.beta,
            "A": [],  # Placeholder A matrix
            "gamma": [],
            "a": self.a,
            "b": self.b,
            "D": self.D,
            "C": self.C,
            "alpha": self.alpha,
            "I": self.I,
            "J": self.J,
            "K": self.K,
            "n": self.n,
            "p": self.p,
            "q": self.q,
            "r": 0,
            "lb": self.lb,
            "ub": self.ub,
            "n_lower_level_constraints": self.n_lower_level_constraints,
            "instance_type": instance_type,
            "instance_file": instance_file,
            "bilinearities": True
        }
        return problem_data

    def add_interdiction_data(self):
        """Add additional data for interdiction instances."""
        # Add interdiction constraints
        for i in range(len(self.D[0])):
            identity_row = [0.0] * len(self.D[0])
            identity_row[i] = 1.0
            self.D.append(identity_row)

        # Add interdiction constraints to C matrix
        original_C = []
        for j in range(self.n_lower_level_constraints):
            original_C.append(
                [0.0] * len(self.D[0])
            )  # Assuming C should be zeros for the first row

        # Add identity matrix rows below
        self.C = original_C  # First rows
        for i in range(len(self.C[0])):
            identity_row = [0.0] * len(self.C[0])
            identity_row[i] = 1.0
            self.C.append(identity_row)

        # Add alpha values for interdiction constraints
        for i in range(len(self.C[0])):
            self.alpha.append(1.0)
