# Branch-and-Cut for Mixed-Integer Linear Decision-Dependent Robust Optimization

## Description

This repository contains the open-source implementations accompanying the paper "Branch-and-Cut for Mixed-Integer Linear Decision-Dependent Robust Optimization" by [Henri Lefebvre](https://henrilefebvre.com/), [Martin Schmidt](https://martinschmidt.squarespace.com/), [Simon Stevens](https://simstevens.github.io/) and [Johannes Thürauf](https://www.johannesthuerauf.com/).

## Prerequisites

The methods are implemented in `Python 3.12.2` and use `CPLEX 22.1.1` as underlying MILP solver. Visit [CPLEX's official website](https://www.ibm.com/de-de/products/ilog-cplex-optimization-studio) for details on how to obtain a license. 
Moreover, all instances can also be solved using `MibS` or `Yasol`. For more details on how to install `MibS` visit the [MibS Quick Start Guide](https://coin-or.github.io/MibS/). For more details on `Yasol` visit the [Yasol Website](https://yasolqipsolver.github.io/yasol.github.io/About_Yasol/) or the [Yasol GitHub](https://github.com/MichaelHartisch/Yasol). The `.mps`, `.aux` and `.qlp` files needed for solving the instances with those solvers are also provided in the `instances` folder.

## Quick start

You can run the branch-and-cut module using the provided script `run_bnc.sh`:
```bash
chmod +x run_bnc.sh           # one-time
./run_bnc.sh --parameter1 value1 --parameter2 value2 ...  # run with parameters
```

## Command-line parameters

All parameters defined in `src/run_bnc.py` are forwarded by `run_bnc.sh`. The parameters without default values are required. The following table summarizes the parameters:

| Flag | Type | Choices | Default | Description |
| --- | --- | --- | --- | --- |
| `--instance_file` | str | – | – | Path to the instance. E.g., `instances/knapsack/knapsack_general/knapsack_100_1_1.kp` |
| `--instance_type` | str | `knapsack`, `bobilib` | – | Instance family. |
| `--cuts` | str | `branchandbound`, `intersection`, `interdiction`, `nogood` | – | Cut strategy. |
| `--lower_level` | str | `general`, `interdiction` | – | Lower-level problem setting. For BOBILib instances, only `general` is supported. |
| `--projected` | int | `0`, `1` | `1` | Use projected formulation (1) or unprojected (0). Projected corresponds to the DDRO framework, unprojected corresponds to the bilevel framework. |
| `--separation` | str | `integer`, `fractional` | `integer` | Separation strategy for intersection cuts. |
| `--time_lim` | float | – | `60.0` | Total time limit in seconds, including the bound-computation time recorded in a colocated conversion manifest. |
| `--output_csv` | str | – | `results/summary.csv` | CSV file to which one summary row is appended after each run. |
| `--verbose_level` | int | `0`, `1`, `2`, `3` | `0` | Verbosity: 0=silent, 1=info, 2=time logging, 3=debug. |
| `--max_cuts` | int | – | `20` | Max cuts per node. |
| `--cplex_cuts` | int | `-1`, `0` | `-1` | Enable CPLEX cuts (0) or disable (-1). |
| `--only_root_node` | int | `0`, `1` | `0` | Separate cuts only at root (1) or at all nodes (0). |
| `--tolerance` | float | – | `1e-6` | Optimality gap tolerance. |
| `--write_lps` | int | `0`, `1` | `0` | Write LP files of node problems and subproblems during solving. |

Notes:
- For BOBILib instances, ensure the instance appears in `instances/bobilib/solvable-instances.csv`. If not, there exists no known solution to the instance.
- If the instance directory contains `conversion_manifest.csv`, the runner reads the matching `bound_time_seconds`, subtracts it from `--time_lim`, and reports both the model and total times in the result CSV. If no manifest is present, the bound time is zero.
