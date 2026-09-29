"""Parser for the original decision-dependent robust TSP instances."""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal
from pathlib import Path


@dataclass(frozen=True)
class SPEdge:
    index: int
    first: int
    second: int
    nominal_cost: Decimal
    score: Decimal
    adjacent_edges: tuple[int, ...]


@dataclass(frozen=True)
class SPInstance:
    path: Path
    name: str
    num_vertices: int
    num_edges: int
    a_0: int
    b_0: tuple[int, ...]
    a: tuple[Decimal, ...]
    b: tuple[dict[int, Decimal], ...]
    b_r_1_lb: tuple[Decimal, ...]
    b_r_1_ub: tuple[Decimal, ...]
    b_r_0_lb: tuple[Decimal, ...]
    b_r_0_ub: tuple[Decimal, ...]
    c_r_1_ub: tuple[Decimal, ...]
    b_r_1_lb_time: Decimal
    b_r_1_ub_time: Decimal
    b_r_0_lb_time: Decimal
    b_r_0_ub_time: Decimal
    edges: tuple[SPEdge, ...]


def _unique_tokens(lines: list[str], key: str) -> list[str]:
    matches = [line.split()[1:] for line in lines if line.split()[:1] == [key]]
    if len(matches) != 1:
        raise ValueError(f"Expected exactly one {key!r} row, found {len(matches)}")
    return matches[0]


def _integer(value: str, label: str) -> int:
    decimal = Decimal(value)
    if decimal != decimal.to_integral_value():
        raise ValueError(f"{label} must be integral, got {decimal}")
    return int(decimal)


def _vector(lines: list[str], key: str, size: int) -> tuple[Decimal, ...]:
    values = tuple(Decimal(value) for value in _unique_tokens(lines, key))
    if len(values) != size:
        raise ValueError(f"{key} contains {len(values)} values; expected {size}")
    return values


def _integer_vector(lines: list[str], key: str, size: int) -> tuple[int, ...]:
    values = tuple(_integer(value, key) for value in _unique_tokens(lines, key))
    if len(values) != size:
        raise ValueError(f"{key} contains {len(values)} values; expected {size}")
    return values


def _section_rows(
    lines: list[str], section: str, next_section: str, expected: int
) -> list[str]:
    try:
        start = lines.index(section) + 1
        end = lines.index(next_section, start)
    except ValueError as exc:
        raise ValueError(f"Missing {section}/{next_section} section boundary") from exc
    rows = lines[start:end]
    if len(rows) != expected:
        raise ValueError(f"{section} contains {len(rows)} rows; expected {expected}")
    return rows


def parse_sp_instance(path: str | Path) -> SPInstance:
    """Parse and validate one original ``InstanceSP`` text file."""
    instance_path = Path(path).resolve()
    lines = [
        line.strip()
        for line in instance_path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]

    name = " ".join(_unique_tokens(lines, "NAME"))
    num_vertices = _integer(_unique_tokens(lines, "NUM_VERTICES")[0], "NUM_VERTICES")
    num_edges = _integer(_unique_tokens(lines, "NUM_EDGES")[0], "NUM_EDGES")
    if num_vertices < 3 or num_edges < num_vertices:
        raise ValueError(
            f"Invalid TSP dimensions: {num_vertices} vertices, {num_edges} edges"
        )

    a_0 = _integer(_unique_tokens(lines, "a_0")[0], "a_0")
    b_0 = _integer_vector(lines, "b_0r", num_edges)
    a = _vector(lines, "a_r", num_edges)
    b_r_1_lb = _vector(lines, "B_r_1_lb", num_edges)
    b_r_1_ub = _vector(lines, "B_r_1_ub", num_edges)
    b_r_0_lb = _vector(lines, "B_r_0_lb", num_edges)
    b_r_0_ub = _vector(lines, "B_r_0_ub", num_edges)
    c_r_1_ub = _vector(lines, "c_r_1_ub", num_edges)
    b_r_1_lb_time = Decimal(_unique_tokens(lines, "B_r_1_lb_time")[0])
    b_r_1_ub_time = Decimal(_unique_tokens(lines, "B_r_1_ub_time")[0])
    b_r_0_lb_time = Decimal(_unique_tokens(lines, "B_r_0_lb_time")[0])
    b_r_0_ub_time = Decimal(_unique_tokens(lines, "B_r_0_ub_time")[0])

    matrix_rows = _section_rows(lines, "b_rm", "EDGES", num_edges)
    sparse_matrix: list[dict[int, Decimal]] = [{} for _ in range(num_edges)]
    for expected_row, row in enumerate(matrix_rows):
        parts = row.split()
        if len(parts) < 2:
            raise ValueError(f"Malformed b_rm row {expected_row}")
        row_id = _integer(parts[0], "b_rm row id")
        nnz = _integer(parts[1], f"b_rm[{row_id}] nnz")
        if row_id != expected_row or len(parts) != nnz + 2:
            raise ValueError(f"Malformed b_rm row {expected_row}")
        for token in parts[2:]:
            try:
                column_token, value_token = token.split(":", 1)
            except ValueError as exc:
                raise ValueError(f"Invalid b_rm token {token!r}") from exc
            column = _integer(column_token, "b_rm column")
            if not 0 <= column < num_edges or column == row_id:
                raise ValueError(f"Invalid b_rm column {column} in row {row_id}")
            if column in sparse_matrix[row_id]:
                raise ValueError(f"Duplicate b_rm entry ({row_id}, {column})")
            value = Decimal(value_token)
            if value < 0:
                raise ValueError("The converter requires nonnegative b_rm coefficients")
            if value:
                sparse_matrix[row_id][column] = value

    edge_rows = _section_rows(lines, "EDGES", "VERTICES", num_edges)
    parsed_edges: list[SPEdge | None] = [None] * num_edges
    undirected_edges: set[tuple[int, int]] = set()
    for expected_edge, row in enumerate(edge_rows):
        parts = row.split()
        if len(parts) < 6:
            raise ValueError(f"Malformed EDGES row {expected_edge}")
        edge_id = _integer(parts[0], "edge id")
        first = _integer(parts[1], f"edge {edge_id} first endpoint")
        second = _integer(parts[2], f"edge {edge_id} second endpoint")
        adjacent_count = _integer(parts[5], f"edge {edge_id} adjacent count")
        adjacent = tuple(_integer(value, "adjacent edge") for value in parts[6:])
        if edge_id != expected_edge or len(adjacent) != adjacent_count:
            raise ValueError(f"Malformed EDGES row {expected_edge}")
        if not (0 <= first < num_vertices and 0 <= second < num_vertices):
            raise ValueError(f"Invalid endpoints for edge {edge_id}: ({first}, {second})")
        if first == second:
            raise ValueError(f"Self-loop at edge {edge_id}")
        key = tuple(sorted((first, second)))
        if key in undirected_edges:
            raise ValueError(f"Duplicate undirected edge {key}")
        undirected_edges.add(key)
        if any(not 0 <= other < num_edges or other == edge_id for other in adjacent):
            raise ValueError(f"Invalid adjacency list for edge {edge_id}")
        parsed_edges[edge_id] = SPEdge(
            edge_id,
            first,
            second,
            Decimal(parts[3]),
            Decimal(parts[4]),
            adjacent,
        )

    if any(value < 0 for value in a):
        raise ValueError("The converter requires nonnegative a_r coefficients")
    bound_vectors = {
        "B_r_1_lb": b_r_1_lb,
        "B_r_1_ub": b_r_1_ub,
        "B_r_0_lb": b_r_0_lb,
        "B_r_0_ub": b_r_0_ub,
    }
    for label, values in bound_vectors.items():
        if any(value < 0 for value in values):
            raise ValueError(f"{label} must be nonnegative")
    for i in range(num_edges):
        if b_r_1_lb[i] > b_r_1_ub[i]:
            raise ValueError(f"Invalid conditional-one bounds for edge {i}")
        if b_r_0_lb[i] > b_r_0_ub[i]:
            raise ValueError(f"Invalid conditional-zero bounds for edge {i}")
    if any(value < 0 for value in c_r_1_ub):
        raise ValueError("c_r_1_ub must be nonnegative")
    bound_times = (
        b_r_1_lb_time,
        b_r_1_ub_time,
        b_r_0_lb_time,
        b_r_0_ub_time,
    )
    if any(value < 0 for value in bound_times):
        raise ValueError("Bound-computation times must be nonnegative")

    return SPInstance(
        path=instance_path,
        name=name,
        num_vertices=num_vertices,
        num_edges=num_edges,
        a_0=a_0,
        b_0=b_0,
        a=a,
        b=tuple(sparse_matrix),
        b_r_1_lb=b_r_1_lb,
        b_r_1_ub=b_r_1_ub,
        b_r_0_lb=b_r_0_lb,
        b_r_0_ub=b_r_0_ub,
        c_r_1_ub=c_r_1_ub,
        b_r_1_lb_time=b_r_1_lb_time,
        b_r_1_ub_time=b_r_1_ub_time,
        b_r_0_lb_time=b_r_0_lb_time,
        b_r_0_ub_time=b_r_0_ub_time,
        edges=tuple(edge for edge in parsed_edges if edge is not None),
    )
