#!/usr/bin/env python3
"""Convert original DDRO TSP text instances to BOBILib MPS/AUX pairs.

Degree-two constraints are written to the MPS.  Connectivity constraints are
intentionally omitted because ``run_bnc_tsp.py`` separates them lazily.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

from convert_kp_to_bobilib import (
    LinearMPSModel,
    render_aux,
    scaled_integer,
    write_text_atomically,
)
from parse_sp import SPInstance, parse_sp_instance


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT_ROOT = REPOSITORY_ROOT / "instances" / "originalSP"
DEFAULT_OUTPUT_DIR = REPOSITORY_ROOT / "instances" / "sp_aggregated_bobilib"
DEFAULT_SCALE = 100
PAPER_SIZES = (25, 50)
PAPER_INSTANCE_COUNT = 40


def instance_identity(instance: SPInstance) -> tuple[int, int]:
    match = re.search(r"_N(\d+).*_ID(\d+)$", instance.path.stem)
    if match is None:
        raise ValueError(f"Cannot extract N/ID from {instance.path.name!r}")
    size, replicate = map(int, match.groups())
    if size != instance.num_vertices:
        raise ValueError(
            f"Filename size {size} differs from NUM_VERTICES {instance.num_vertices}"
        )
    return size, replicate


def build_bilevel_model(
    instance: SPInstance, scale: int
) -> tuple[LinearMPSModel, list[str], list[str], int]:
    n = instance.num_edges
    model = LinearMPSModel(instance.path.stem)
    x_names = [f"x_{i}" for i in range(n)]
    z_names = [f"z_{i}" for i in range(n)]
    r_names = [f"r_{i}" for i in range(n)]
    s_names = [f"s_{i}" for i in range(n)]
    u_names = [f"u_{i}" for i in range(n)]
    p_names = [f"p_{i}" for i in range(n)]

    nominal_costs = [
        scaled_integer(edge.nominal_cost, scale, f"nominal_cost[{edge.index}]")
        for edge in instance.edges
    ]
    scaled_a = [
        scaled_integer(value, scale, f"a_r[{i}]")
        for i, value in enumerate(instance.a)
    ]
    scaled_b = [
        {
            j: scaled_integer(value, scale, f"b_rm[{i},{j}]")
            for j, value in row.items()
        }
        for i, row in enumerate(instance.b)
    ]
    scaled_b_r_1_lb = [
        scaled_integer(value, scale, f"B_r_1_lb[{i}]")
        for i, value in enumerate(instance.b_r_1_lb)
    ]
    scaled_b_r_1_ub = [
        scaled_integer(value, scale, f"B_r_1_ub[{i}]")
        for i, value in enumerate(instance.b_r_1_ub)
    ]
    scaled_b_r_0_lb = [
        scaled_integer(value, scale, f"B_r_0_lb[{i}]")
        for i, value in enumerate(instance.b_r_0_lb)
    ]
    scaled_b_r_0_ub = [
        scaled_integer(value, scale, f"B_r_0_ub[{i}]")
        for i, value in enumerate(instance.b_r_0_ub)
    ]
    scaled_z_upper_bounds = [
        scaled_integer(value, scale, f"c_r_1_ub[{i}]")
        for i, value in enumerate(instance.c_r_1_ub)
    ]

    for name in x_names:
        model.add_variable(name, binary=True, upper=1)
    for i, name in enumerate(z_names):
        model.add_variable(name, binary=False, lower=0, upper=scaled_z_upper_bounds[i])
    for i, name in enumerate(r_names):
        model.add_variable(name, binary=False, lower=0, upper=scaled_b_r_1_ub[i])
    for i, name in enumerate(s_names):
        lower = min(scaled_b_r_1_lb[i], scaled_b_r_0_lb[i])
        upper = max(scaled_b_r_1_ub[i], scaled_b_r_0_ub[i])
        model.add_variable(name, binary=False, lower=lower, upper=upper)
    for name in u_names:
        model.add_variable(name, binary=True, upper=1)
    for i, name in enumerate(p_names):
        model.add_variable(name, binary=False, lower=0, upper=scaled_z_upper_bounds[i])

    for i, name in enumerate(x_names):
        model.add_objective_coefficient(name, nominal_costs[i])
    for name in p_names:
        model.add_objective_coefficient(name, 1)

    for vertex in range(instance.num_vertices):
        incident = (
            (x_names[edge.index], 1)
            for edge in instance.edges
            if edge.first == vertex or edge.second == vertex
        )
        model.add_row(f"UL_degree_{vertex}", "E", 2, incident)

    # Paper-consistent aggregated linearization:
    #   s_i = sum_{j != i} b_ij x_j,
    #   r_i = x_i s_i,
    #   z_i = a_i x_i + r_i.
    # All four stored conditional bounds are retained so that the initial LP
    # relaxation matches the aggregated TSP formulation used in the paper.
    for i in range(n):
        model.add_row(
            f"UL_s_{i}",
            "E",
            0,
            [
                (s_names[i], 1),
                *((x_names[j], -value) for j, value in scaled_b[i].items()),
            ],
        )
        model.add_row(
            f"UL_z_{i}",
            "E",
            0,
            (
                (z_names[i], 1),
                (x_names[i], -scaled_a[i]),
                (r_names[i], -1),
            ),
        )
        model.add_row(
            f"UL_r1lo_{i}",
            "G",
            0,
            ((r_names[i], 1), (x_names[i], -scaled_b_r_1_lb[i])),
        )
        model.add_row(
            f"UL_r1up_{i}",
            "L",
            0,
            ((r_names[i], 1), (x_names[i], -scaled_b_r_1_ub[i])),
        )
        model.add_row(
            f"UL_r0lo_{i}",
            "G",
            -scaled_b_r_0_ub[i],
            (
                (r_names[i], 1),
                (s_names[i], -1),
                (x_names[i], -scaled_b_r_0_ub[i]),
            ),
        )
        model.add_row(
            f"UL_r0up_{i}",
            "L",
            -scaled_b_r_0_lb[i],
            (
                (r_names[i], 1),
                (s_names[i], -1),
                (x_names[i], -scaled_b_r_0_lb[i]),
            ),
        )

    lower_rows: list[str] = []
    budget_terms = [(name, 1) for name in u_names]
    budget_terms.extend((x_names[i], -instance.b_0[i]) for i in range(n))
    model.add_row("LL_budget", "L", instance.a_0, budget_terms)
    lower_rows.append("LL_budget")

    for i in range(n):
        upper = scaled_z_upper_bounds[i]
        rows = (
            (f"LL_pub_{i}", "L", 0, ((p_names[i], 1), (u_names[i], -upper))),
            (f"LL_pz_{i}", "L", 0, ((p_names[i], 1), (z_names[i], -1))),
            (
                f"LL_plo_{i}",
                "G",
                -upper,
                ((p_names[i], 1), (z_names[i], -1), (u_names[i], -upper)),
            ),
        )
        for row_name, sense, rhs, terms in rows:
            model.add_row(row_name, sense, rhs, terms)
            lower_rows.append(row_name)

    nonzero_interactions = sum(len(row) for row in instance.b)
    return model, u_names + p_names, lower_rows, nonzero_interactions


def _relative_or_absolute(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(REPOSITORY_ROOT))
    except ValueError:
        return str(path.resolve())


def discover_instances(input_root: Path, sizes: set[int]) -> list[SPInstance]:
    instances = [parse_sp_instance(path) for path in sorted(input_root.rglob("*.txt"))]
    selected = [item for item in instances if item.num_vertices in sizes]
    selected.sort(key=lambda item: (*instance_identity(item), item.path.name))
    return selected


def convert_instances(args: argparse.Namespace) -> None:
    input_root = args.input_root.resolve()
    output_dir = args.output_dir.resolve()
    if not input_root.is_dir():
        raise FileNotFoundError(f"Input directory does not exist: {input_root}")
    instances = discover_instances(input_root, set(args.sizes))
    if args.expected_count and len(instances) != args.expected_count:
        raise ValueError(f"Selected {len(instances)} SP instances; expected {args.expected_count}")

    print(f"Selected {len(instances)} SP instances from {input_root}")
    if args.dry_run:
        for instance in instances:
            size, replicate = instance_identity(instance)
            interaction_count = sum(len(row) for row in instance.b)
            print(
                f"  {instance.path.name}: vertices={size}, edges={instance.num_edges}, "
                f"interactions={interaction_count}, replicate={replicate}"
            )
        return

    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_rows: list[dict[str, str | int]] = []
    for instance in instances:
        size, replicate = instance_identity(instance)
        model, lower_variables, lower_rows, interaction_count = build_bilevel_model(
            instance, args.scale
        )
        stem = instance.path.stem
        mps_path = output_dir / f"{stem}.mps"
        aux_path = output_dir / f"{stem}.aux"
        write_text_atomically(mps_path, model.render(), args.overwrite)
        write_text_atomically(
            aux_path,
            render_aux(stem, lower_variables, lower_rows),
            args.overwrite,
        )
        manifest_rows.append(
            {
                "instance": stem,
                "source_file": _relative_or_absolute(instance.path),
                "problem": "TSP",
                "formulation": "DDRO_BI",
                "linearcons": "A",
                "n": size,
                "num_edges": instance.num_edges,
                "density": "",
                "replicate": replicate,
                "scale": args.scale,
                "input_objective_sense": "min",
                "reported_objective_sense": "min",
                "objective_scale": args.scale,
                "objective_multiplier": 1.0 / args.scale,
                "bound_time_seconds": str(
                    instance.b_r_1_lb_time
                    + instance.b_r_1_ub_time
                    + instance.b_r_0_lb_time
                    + instance.b_r_0_ub_time
                ),
                "upper_variables": 4 * instance.num_edges,
                "lower_variables": len(lower_variables),
                "upper_constraints": size + 6 * instance.num_edges,
                "lower_constraints": len(lower_rows),
                "nonzero_interactions": interaction_count,
                "mps_file": mps_path.name,
                "aux_file": aux_path.name,
            }
        )

    manifest_path = output_dir / "conversion_manifest.csv"
    if manifest_path.exists() and not args.overwrite:
        raise FileExistsError(f"Refusing to overwrite existing file: {manifest_path}")
    temporary = manifest_path.with_name(f".{manifest_path.name}.tmp")
    with temporary.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(manifest_rows[0]))
        writer.writeheader()
        writer.writerows(manifest_rows)
    temporary.replace(manifest_path)
    print(f"Wrote {len(instances)} MPS/AUX pairs to {output_dir}")
    print(f"Wrote conversion metadata to {manifest_path}")


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert original DDRO TSP text instances to BOBILib MPS/AUX pairs."
    )
    parser.add_argument("--input-root", type=Path, default=DEFAULT_INPUT_ROOT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--sizes", type=int, nargs="+", default=list(PAPER_SIZES))
    parser.add_argument("--scale", type=int, default=DEFAULT_SCALE)
    parser.add_argument("--expected-count", type=int, default=PAPER_INSTANCE_COUNT)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    if args.scale <= 0:
        parser.error("--scale must be positive")
    if args.expected_count < 0:
        parser.error("--expected-count must be nonnegative")
    return args


if __name__ == "__main__":
    convert_instances(parse_arguments())
