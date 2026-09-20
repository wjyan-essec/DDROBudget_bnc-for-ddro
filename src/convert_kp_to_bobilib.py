#!/usr/bin/env python3
"""Convert the project's disaggregated KP instances to BOBILib MPS/AUX pairs.

The generated files encode the original bilevel model (28) from the DDRO
branch-and-cut paper.  The existing Python implementation can then construct
the projected DDRO model (29)--(30) by running with ``--projected 1``.

No optimization package is needed to generate the files.  All decimal profit
and deviation data are scaled by 100 so that the BOBILib input is pure integer,
as required by the current implementation.
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import OrderedDict
from dataclasses import dataclass
from decimal import Decimal
from pathlib import Path
from typing import Iterable


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT_ROOT = REPOSITORY_ROOT / "instances" / "originalKP"
DEFAULT_OUTPUT_DIR = REPOSITORY_ROOT / "instances" / "kp_disaggregated_bobilib"
DEFAULT_SCALE = 100
PAPER_SIZES = (50, 100)
PAPER_DENSITIES = (50, 75)
PAPER_INSTANCE_COUNT = 40


@dataclass(frozen=True)
class KPInstance:
    path: Path
    name: str
    source: str
    size: int
    density: int
    replicate: int
    nominal_profits: tuple[Decimal, ...]
    weights: tuple[int, ...]
    capacity: int
    a_0: int
    b_0: tuple[int, ...]
    a: tuple[Decimal, ...]
    b: tuple[tuple[Decimal, ...], ...]
    c_i_1_ub: tuple[Decimal, ...]
    b_i_1_ub_time: Decimal


@dataclass(frozen=True)
class VariableBound:
    name: str
    binary: bool
    lower: int = 0
    upper: int = 0


class LinearMPSModel:
    """Small deterministic writer for the subset of free MPS used here."""

    def __init__(self, name: str) -> None:
        self.name = name
        self.rows: list[tuple[str, str, int]] = []
        self.row_names: set[str] = set()
        self.columns: OrderedDict[str, OrderedDict[str, int]] = OrderedDict()
        self.bounds: OrderedDict[str, VariableBound] = OrderedDict()

    def add_variable(
        self, name: str, *, binary: bool, lower: int = 0, upper: int = 0
    ) -> None:
        if name in self.columns:
            raise ValueError(f"Duplicate variable name: {name}")
        if upper < lower:
            raise ValueError(f"Invalid bounds for {name}: [{lower}, {upper}]")
        self.columns[name] = OrderedDict()
        self.bounds[name] = VariableBound(name, binary, lower, upper)

    def add_objective_coefficient(self, variable: str, coefficient: int) -> None:
        self._add_coefficient(variable, "OBJ", coefficient)

    def add_row(
        self,
        name: str,
        sense: str,
        rhs: int,
        terms: Iterable[tuple[str, int]],
    ) -> None:
        if sense not in {"L", "E", "G"}:
            raise ValueError(f"Unsupported row sense {sense!r} for {name}")
        if name in self.row_names or name == "OBJ":
            raise ValueError(f"Duplicate row name: {name}")
        self.rows.append((sense, name, rhs))
        self.row_names.add(name)
        for variable, coefficient in terms:
            self._add_coefficient(variable, name, coefficient)

    def _add_coefficient(self, variable: str, row: str, coefficient: int) -> None:
        if variable not in self.columns:
            raise ValueError(f"Unknown variable {variable!r} in row {row!r}")
        if coefficient == 0:
            return
        existing = self.columns[variable].get(row, 0)
        updated = existing + coefficient
        if updated:
            self.columns[variable][row] = updated
        elif row in self.columns[variable]:
            del self.columns[variable][row]

    def render(self) -> str:
        lines = [f"NAME          {self.name}", "OBJSENSE", " MIN", "ROWS", " N  OBJ"]
        lines.extend(f" {sense}  {name}" for sense, name, _ in self.rows)
        lines.append("COLUMNS")
        lines.append("    MARK0000  'MARKER'                 'INTORG'")
        for variable, entries in self.columns.items():
            for row, coefficient in entries.items():
                lines.append(f"    {variable}  {row}  {coefficient}")
        lines.append("    MARK0001  'MARKER'                 'INTEND'")
        lines.append("RHS")
        for _, row, rhs in self.rows:
            if rhs:
                lines.append(f"    RHS1  {row}  {rhs}")
        lines.append("BOUNDS")
        for bound in self.bounds.values():
            if bound.binary:
                lines.append(f" BV BND1  {bound.name}")
            else:
                lines.append(f" LO BND1  {bound.name}  {bound.lower}")
                lines.append(f" UP BND1  {bound.name}  {bound.upper}")
        lines.append("ENDATA")
        return "\n".join(lines) + "\n"


def _tokens_for_key(lines: list[str], key: str) -> list[str]:
    matches = [line.split()[1:] for line in lines if line.split() and line.split()[0] == key]
    if len(matches) != 1:
        raise ValueError(f"Expected exactly one {key!r} row, found {len(matches)}")
    return matches[0]


def _decimal_vector(lines: list[str], key: str, size: int) -> tuple[Decimal, ...]:
    values = tuple(Decimal(value) for value in _tokens_for_key(lines, key))
    if len(values) != size:
        raise ValueError(f"{key} contains {len(values)} values; expected {size}")
    return values


def _integer(value: str | Decimal, label: str) -> int:
    decimal_value = value if isinstance(value, Decimal) else Decimal(value)
    if decimal_value != decimal_value.to_integral_value():
        raise ValueError(f"{label} must be integral, got {decimal_value}")
    return int(decimal_value)


def _integer_vector(lines: list[str], key: str, size: int) -> tuple[int, ...]:
    values = tuple(_integer(value, key) for value in _tokens_for_key(lines, key))
    if len(values) != size:
        raise ValueError(f"{key} contains {len(values)} values; expected {size}")
    return values


def parse_kp_instance(path: Path) -> KPInstance:
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    name = " ".join(_tokens_for_key(lines, "NAME"))
    size = _integer(_tokens_for_key(lines, "NUM_ITEMS")[0], "NUM_ITEMS")
    source = _tokens_for_key(lines, "SOURCE")[0]
    source_match = re.fullmatch(r"(?:jeu|r)_(\d+)_(\d+)_(\d+)", source)
    if source_match is None:
        raise ValueError(f"Cannot extract size/density/replicate from SOURCE {source!r}")
    source_size, density, replicate = map(int, source_match.groups())
    if source_size != size:
        raise ValueError(f"SOURCE size {source_size} differs from NUM_ITEMS {size}")

    nominal_profits = _decimal_vector(lines, "nominal_profits", size)
    weights = _integer_vector(lines, "weights", size)
    capacity = _integer(_tokens_for_key(lines, "capacity")[0], "capacity")
    a_0 = _integer(_tokens_for_key(lines, "a_0")[0], "a_0")
    b_0 = _integer_vector(lines, "b_0i", size)
    a = _decimal_vector(lines, "a_i", size)
    c_i_1_ub = _decimal_vector(lines, "c_i_1_ub", size)
    b_i_1_ub_time = Decimal(_tokens_for_key(lines, "B_i_1_ub_time")[0])

    try:
        matrix_start = lines.index("b_ij") + 1
    except ValueError as exc:
        raise ValueError("Missing b_ij matrix") from exc
    matrix: list[tuple[Decimal, ...]] = []
    for expected_row, line in enumerate(lines[matrix_start : matrix_start + size]):
        parts = line.split()
        if len(parts) != size + 1:
            raise ValueError(
                f"b_ij row {expected_row} has {len(parts) - 1} entries; expected {size}"
            )
        row_id = _integer(parts[0], f"b_ij row id {expected_row}")
        if row_id != expected_row:
            raise ValueError(f"Expected b_ij row {expected_row}, found row {row_id}")
        matrix.append(tuple(Decimal(value) for value in parts[1:]))
    if len(matrix) != size:
        raise ValueError(f"b_ij contains {len(matrix)} rows; expected {size}")

    if any(value < 0 for value in a):
        raise ValueError("The converter requires nonnegative a_i coefficients")
    if any(value < 0 for row in matrix for value in row):
        raise ValueError("The converter requires nonnegative b_ij coefficients")
    if any(matrix[i][i] != 0 for i in range(size)):
        raise ValueError("The converter requires b_ii = 0")

    return KPInstance(
        path=path,
        name=name,
        source=source,
        size=size,
        density=density,
        replicate=replicate,
        nominal_profits=nominal_profits,
        weights=weights,
        capacity=capacity,
        a_0=a_0,
        b_0=b_0,
        a=a,
        b=tuple(matrix),
        c_i_1_ub=c_i_1_ub,
        b_i_1_ub_time=b_i_1_ub_time,
    )


def scaled_integer(value: Decimal, scale: int, label: str) -> int:
    scaled = value * scale
    if scaled != scaled.to_integral_value():
        raise ValueError(
            f"{label}={value} is not integral after multiplication by scale {scale}"
        )
    return int(scaled)


def build_bilevel_model(
    instance: KPInstance, scale: int
) -> tuple[LinearMPSModel, list[str], list[str]]:
    n = instance.size
    model_name = instance.path.stem
    model = LinearMPSModel(model_name)

    x_names = [f"x_{i}" for i in range(n)]
    z_names = [f"z_{i}" for i in range(n)]
    pair_names = {(i, j): f"t_{i}_{j}" for i in range(n) for j in range(i + 1, n)}
    u_names = [f"u_{i}" for i in range(n)]
    p_names = [f"p_{i}" for i in range(n)]

    scaled_profits = [
        scaled_integer(value, scale, f"nominal_profits[{i}]")
        for i, value in enumerate(instance.nominal_profits)
    ]
    scaled_a = [
        scaled_integer(value, scale, f"a_i[{i}]") for i, value in enumerate(instance.a)
    ]
    scaled_b = [
        [scaled_integer(value, scale, f"b_ij[{i},{j}]") for j, value in enumerate(row)]
        for i, row in enumerate(instance.b)
    ]
    scaled_upper_bounds = [
        scaled_integer(value, scale, f"c_i_1_ub[{i}]")
        for i, value in enumerate(instance.c_i_1_ub)
    ]

    for name in x_names:
        model.add_variable(name, binary=True, upper=1)
    for i, name in enumerate(z_names):
        model.add_variable(name, binary=False, lower=0, upper=scaled_upper_bounds[i])
    for name in pair_names.values():
        model.add_variable(name, binary=True, upper=1)
    for name in u_names:
        model.add_variable(name, binary=True, upper=1)
    for i, name in enumerate(p_names):
        model.add_variable(name, binary=False, lower=0, upper=scaled_upper_bounds[i])

    for i, name in enumerate(x_names):
        model.add_objective_coefficient(name, -scaled_profits[i])
    for name in p_names:
        model.add_objective_coefficient(name, 1)

    model.add_row(
        "UL_capacity",
        "L",
        instance.capacity,
        ((x_names[i], instance.weights[i]) for i in range(n)),
    )

    for i in range(n):
        z_terms: list[tuple[str, int]] = [(z_names[i], 1), (x_names[i], -scaled_a[i])]
        for j in range(n):
            if i == j:
                continue
            pair = (j, i) if j < i else (i, j)
            z_terms.append((pair_names[pair], -scaled_b[i][j]))
        model.add_row(f"UL_z_{i}", "E", 0, z_terms)

    for (i, j), t_name in pair_names.items():
        model.add_row(f"UL_tux_{i}_{j}", "L", 0, ((t_name, 1), (x_names[i], -1)))
        model.add_row(f"UL_tuy_{i}_{j}", "L", 0, ((t_name, 1), (x_names[j], -1)))
        model.add_row(
            f"UL_tlo_{i}_{j}",
            "G",
            -1,
            ((t_name, 1), (x_names[i], -1), (x_names[j], -1)),
        )

    lower_rows: list[str] = []
    budget_row = "LL_budget"
    budget_terms = [(name, 1) for name in u_names]
    budget_terms.extend((x_names[i], -instance.b_0[i]) for i in range(n))
    model.add_row(budget_row, "L", instance.a_0, budget_terms)
    lower_rows.append(budget_row)

    for i in range(n):
        upper = scaled_upper_bounds[i]
        product_rows = (
            (f"LL_pub_{i}", "L", 0, ((p_names[i], 1), (u_names[i], -upper))),
            (f"LL_pz_{i}", "L", 0, ((p_names[i], 1), (z_names[i], -1))),
            (
                f"LL_plo_{i}",
                "G",
                -upper,
                ((p_names[i], 1), (z_names[i], -1), (u_names[i], -upper)),
            ),
        )
        for row_name, sense, rhs, terms in product_rows:
            model.add_row(row_name, sense, rhs, terms)
            lower_rows.append(row_name)

    lower_variables = u_names + p_names
    if len(lower_variables) != 2 * n or len(lower_rows) != 3 * n + 1:
        raise AssertionError("Unexpected lower-level dimensions")
    return model, lower_variables, lower_rows


def render_aux(
    model_name: str, lower_variables: list[str], lower_rows: list[str]
) -> str:
    lines = [
        "@NUMVARS",
        str(len(lower_variables)),
        "@NUMCONSTRS",
        str(len(lower_rows)),
        "@VARSBEGIN",
    ]
    for variable in lower_variables:
        coefficient = -1 if variable.startswith("p_") else 0
        lines.append(f"{variable}  {coefficient}")
    lines.extend(["@VARSEND", "@CONSTRSBEGIN"])
    lines.extend(lower_rows)
    lines.extend(
        [
            "@CONSTRSEND",
            "@NAME",
            model_name,
            "@MPS",
            f"{model_name}.mps",
        ]
    )
    return "\n".join(lines) + "\n"


def write_text_atomically(path: Path, content: str, overwrite: bool) -> None:
    if path.exists() and not overwrite:
        raise FileExistsError(f"Refusing to overwrite existing file: {path}")
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(content, encoding="utf-8")
    temporary.replace(path)


def discover_instances(
    input_root: Path, sizes: set[int], densities: set[int]
) -> tuple[list[KPInstance], list[Path]]:
    selected: list[KPInstance] = []
    ignored: list[Path] = []
    for path in sorted(input_root.rglob("*.txt")):
        instance = parse_kp_instance(path)
        if instance.size in sizes and instance.density in densities:
            selected.append(instance)
        else:
            ignored.append(path)
    selected.sort(key=lambda item: (item.size, item.density, item.replicate, item.path.name))
    return selected, ignored


def convert_instances(args: argparse.Namespace) -> None:
    input_root = args.input_root.resolve()
    output_dir = args.output_dir.resolve()
    if not input_root.is_dir():
        raise FileNotFoundError(f"Input directory does not exist: {input_root}")

    instances, ignored = discover_instances(input_root, set(args.sizes), set(args.densities))
    if args.expected_count and len(instances) != args.expected_count:
        raise ValueError(
            f"Selected {len(instances)} instances, expected {args.expected_count}. "
            f"Sizes={args.sizes}, densities={args.densities}."
        )

    print(f"Selected {len(instances)} KP instances from {input_root}")
    print(f"Ignored {len(ignored)} instances outside the requested sizes/densities")
    if args.dry_run:
        for instance in instances:
            print(
                f"  {instance.path.name}: n={instance.size}, "
                f"density={instance.density}, replicate={instance.replicate}"
            )
        return

    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_rows: list[dict[str, str | int]] = []
    for instance in instances:
        model, lower_variables, lower_rows = build_bilevel_model(instance, args.scale)
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
                "source_file": str(instance.path.relative_to(REPOSITORY_ROOT)),
                "n": instance.size,
                "density": instance.density,
                "replicate": instance.replicate,
                "scale": args.scale,
                "bound_time_seconds": str(instance.b_i_1_ub_time),
                "upper_variables": instance.size * 2 + instance.size * (instance.size - 1) // 2,
                "lower_variables": len(lower_variables),
                "upper_constraints": 1
                + instance.size
                + 3 * instance.size * (instance.size - 1) // 2,
                "lower_constraints": len(lower_rows),
                "mps_file": mps_path.name,
                "aux_file": aux_path.name,
            }
        )

    manifest_path = output_dir / "conversion_manifest.csv"
    if manifest_path.exists() and not args.overwrite:
        raise FileExistsError(f"Refusing to overwrite existing file: {manifest_path}")
    temporary_manifest = manifest_path.with_name(f".{manifest_path.name}.tmp")
    with temporary_manifest.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(manifest_rows[0]))
        writer.writeheader()
        writer.writerows(manifest_rows)
    temporary_manifest.replace(manifest_path)

    print(f"Wrote {len(instances)} MPS/AUX pairs to {output_dir}")
    print(f"Wrote conversion metadata to {manifest_path}")


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert disaggregated KP text instances to BOBILib MPS/AUX pairs."
    )
    parser.add_argument("--input-root", type=Path, default=DEFAULT_INPUT_ROOT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--sizes", type=int, nargs="+", default=list(PAPER_SIZES))
    parser.add_argument("--densities", type=int, nargs="+", default=list(PAPER_DENSITIES))
    parser.add_argument("--scale", type=int, default=DEFAULT_SCALE)
    parser.add_argument(
        "--expected-count",
        type=int,
        default=PAPER_INSTANCE_COUNT,
        help="Fail unless this many instances are selected; use 0 to disable the check.",
    )
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
