"""Pure graph helpers used by the TSP connectivity-cut callback."""

from dataclasses import dataclass


@dataclass(frozen=True)
class TSPEdge:
    variable: str
    first: int
    second: int


@dataclass(frozen=True)
class TSPTopology:
    num_vertices: int
    edges: tuple[TSPEdge, ...]

    def validate(self):
        if self.num_vertices < 3:
            raise ValueError("A TSP graph must contain at least three vertices.")
        variables = set()
        undirected_edges = set()
        for edge in self.edges:
            if not 0 <= edge.first < self.num_vertices:
                raise ValueError(f"Invalid endpoint {edge.first} for {edge.variable}.")
            if not 0 <= edge.second < self.num_vertices:
                raise ValueError(f"Invalid endpoint {edge.second} for {edge.variable}.")
            if edge.first == edge.second:
                raise ValueError(f"Self-loop is not allowed: {edge.variable}.")
            if edge.variable in variables:
                raise ValueError(f"Duplicate edge variable: {edge.variable}.")
            variables.add(edge.variable)
            key = tuple(sorted((edge.first, edge.second)))
            if key in undirected_edges:
                raise ValueError(f"Duplicate undirected edge: {key}.")
            undirected_edges.add(key)


def connected_components(num_vertices, edges, edge_values, threshold=0.5):
    """Return components of the support graph induced by selected edges."""
    adjacency = [[] for _ in range(num_vertices)]
    for edge, value in zip(edges, edge_values):
        if value > threshold:
            adjacency[edge.first].append(edge.second)
            adjacency[edge.second].append(edge.first)

    components = []
    seen = [False] * num_vertices
    for start in range(num_vertices):
        if seen[start]:
            continue
        stack = [start]
        seen[start] = True
        component = []
        while stack:
            vertex = stack.pop()
            component.append(vertex)
            for neighbor in adjacency[vertex]:
                if not seen[neighbor]:
                    seen[neighbor] = True
                    stack.append(neighbor)
        components.append(tuple(sorted(component)))
    return components


def canonical_component_shores(num_vertices, components):
    """Return unique shores for component cutsets, identifying complements."""
    all_vertices = frozenset(range(num_vertices))
    shores = set()
    for component in components:
        vertex_set = frozenset(component)
        if not vertex_set or vertex_set == all_vertices:
            continue
        complement = all_vertices - vertex_set
        shore = min(
            (tuple(sorted(vertex_set)), tuple(sorted(complement))),
            key=lambda item: (len(item), item),
        )
        shores.add(shore)
    return sorted(shores, key=lambda item: (len(item), item))


def is_hamiltonian_tour(num_vertices, edges, edge_values, threshold=0.5):
    """Check connectivity and degree two for an undirected selected-edge set."""
    degrees = [0] * num_vertices
    for edge, value in zip(edges, edge_values):
        if value > threshold:
            degrees[edge.first] += 1
            degrees[edge.second] += 1
    return all(degree == 2 for degree in degrees) and len(
        connected_components(num_vertices, edges, edge_values, threshold)
    ) == 1
