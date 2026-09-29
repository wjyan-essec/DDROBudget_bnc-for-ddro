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
- For converted disaggregated KP instances, `incumbent` and `bound` are reported in the original KP maximization scale. The BOBILib minimization values remain available as `raw_incumbent` and `raw_bound`; `objective_scale` records the conversion scale, and `gap_percent` is `100 * gap`.
- The runner checks an existing result CSV header before solving. After a result-schema change, use a new `--output_csv` path (or move the old CSV) instead of appending incompatible rows.

## TSP copy with connectivity cuts

The files under `instances/originalSP` use the companion C++ `InstanceSP` text
format.  Convert them to the integer BOBILib representation used by this
Python framework once:

```bash
python3 src/convert_sp_to_bobilib.py
```

The converter parses the graph, nominal costs, decision-dependent uncertainty
budget, sparse deviation matrix, and all four precomputed conditional-bound
families. Its MPS uses the paper's aggregated `(x,z,r,s)` formulation and
contains the degree-two constraints, but deliberately omits connectivity
constraints. The conversion manifest records the original SP source file and
the sum of all four bound-computation times, so no separate graph file or
edge-variable mapping is needed at solve time.

`run_bnc_tsp.sh` then reuses the MPS/AUX model and DDRO separation code while
adding mandatory TSP connectivity cuts at integer lazy-callback candidates. It
supports two comparable methods:

- `--cuts branchandbound`: `BB-CONN` (connectivity cuts only);
- `--cuts intersection`: `BC-INT-CONN` (connectivity plus intersection cuts).

Connectivity cuts are never subject to `--max_cuts` or `--only_root_node`.
Those options continue to govern intersection separation only.

The runner resolves the original file through the colocated
`conversion_manifest.csv` and maps edge id `e` to upper-level MPS variable
`x_e`. Every mapped edge variable is checked to have bounds in `[0,1]`; a final
incumbent is independently checked for degree two and connectivity.

The TSP callback adds one CPLEX Callable Library symbol to the local Cython
bridge.  Rebuild `src/python_to_cplex_c_api/cplex_c_api_wrapper.pyx` with the
same Python interpreter used to run the experiments before the first TSP run.

```bash
cd src/python_to_cplex_c_api
python setup.py build_ext --inplace
```

Example:

```bash
./run_bnc_tsp.sh \
  --instance_file instances/sp_aggregated_bobilib/SP_N25_rc_1_rGamma_1_dc_0.035_dGamma_0.0005_SP_N25_ID01.mps \
  --cuts branchandbound \
  --time_lim 3600
```

On Windows, the two N25 batch drivers run all 20 converted instances and
write separate logs, `summary.csv`, and a comparison against the C++ EU-A
reference objectives under a fresh `test_runs` directory:

```bat
test_sp_n25_BB_CONN.bat
test_sp_n25_BC_INT_CONN.bat
```

The optional first argument changes the total time limit per instance, for
example `test_sp_n25_BB_CONN.bat 600`. The scripts use
`.venv\Scripts\python.exe` when it exists; otherwise they use `python`, or the
interpreter supplied in `PYTHON_EXE`.
