from __future__ import annotations

import heapq
from collections import deque
from typing import Dict, Iterable, Iterator, List, Optional, Set, Tuple

import numpy as np


class NetworkXNoPath(Exception):
    """Compatibility exception for callers that expect networkx-style errors."""


class _MatrixWrapper:
    def __init__(self, array: np.ndarray):
        self._array = np.array(array, dtype=float, copy=True)

    def toarray(self) -> np.ndarray:
        return np.array(self._array, dtype=float, copy=True)


class Graph:
    def __init__(self):
        self._adj: Dict[int, Dict[int, Dict[str, float]]] = {}

    def __contains__(self, node: int) -> bool:
        return node in self._adj

    def __getitem__(self, node: int) -> Dict[int, Dict[str, float]]:
        return self._adj[node]

    def add_node(self, node: int) -> None:
        self._adj.setdefault(node, {})

    def add_nodes_from(self, nodes: Iterable[int]) -> None:
        for node in nodes:
            self.add_node(int(node))

    def add_edge(self, u: int, v: int, **attrs) -> None:
        iu = int(u)
        iv = int(v)
        self.add_node(iu)
        self.add_node(iv)
        data = dict(self._adj[iu].get(iv, {}))
        data.update(attrs)
        self._adj[iu][iv] = data
        if iu != iv:
            self._adj[iv][iu] = dict(data)

    def has_edge(self, u: int, v: int) -> bool:
        return int(v) in self._adj.get(int(u), {})

    def remove_edge(self, u: int, v: int) -> None:
        iu = int(u)
        iv = int(v)
        self._adj.get(iu, {}).pop(iv, None)
        if iu != iv:
            self._adj.get(iv, {}).pop(iu, None)

    def neighbors(self, node: int) -> Iterator[int]:
        return iter(self._adj.get(int(node), {}))

    def nodes(self) -> List[int]:
        return list(self._adj.keys())

    def edges(self) -> List[Tuple[int, int]]:
        out: List[Tuple[int, int]] = []
        seen: Set[Tuple[int, int]] = set()
        for u, nbrs in self._adj.items():
            for v in nbrs:
                key = (u, v) if u <= v else (v, u)
                if key in seen:
                    continue
                seen.add(key)
                out.append(key)
        return out

    def number_of_nodes(self) -> int:
        return len(self._adj)

    def number_of_edges(self) -> int:
        return len(self.edges())

    def copy(self) -> "Graph":
        out = self.__class__()
        out.add_nodes_from(self.nodes())
        for u, v in self.edges():
            out.add_edge(u, v, **dict(self._adj[u][v]))
        return out

    def subgraph(self, nodes: Iterable[int]) -> "Graph":
        keep = {int(node) for node in nodes}
        out = self.__class__()
        out.add_nodes_from(sorted(keep))
        for u, v in self.edges():
            if u in keep and v in keep:
                out.add_edge(u, v, **dict(self._adj[u][v]))
        return out


class MultiGraph(Graph):
    def __init__(self):
        super().__init__()
        self._edge_mult: Dict[Tuple[int, int], int] = {}

    def add_edge(self, u: int, v: int, **attrs) -> None:
        iu = int(u)
        iv = int(v)
        key = (iu, iv) if iu <= iv else (iv, iu)
        self._edge_mult[key] = self._edge_mult.get(key, 0) + 1
        if self.has_edge(iu, iv):
            new_weight = float(attrs.get("weight", self._adj[iu][iv].get("weight", 1.0)) or 1.0)
            old_weight = float(self._adj[iu][iv].get("weight", new_weight) or new_weight)
            if new_weight < old_weight:
                super().add_edge(iu, iv, **attrs)
            return
        super().add_edge(iu, iv, **attrs)

    def number_of_edges(self) -> int:
        return int(sum(self._edge_mult.values()))

    def copy(self) -> "MultiGraph":
        out = MultiGraph()
        out.add_nodes_from(self.nodes())
        for key, count in self._edge_mult.items():
            u, v = key
            attrs = dict(self._adj[u][v])
            for _ in range(int(count)):
                out.add_edge(u, v, **attrs)
        return out

    def multiplicity(self, u: int, v: int) -> int:
        key = (int(u), int(v))
        if key[0] > key[1]:
            key = (key[1], key[0])
        return int(self._edge_mult.get(key, 0))


def _edge_weight(graph: Graph, u: int, v: int, weight: Optional[str]) -> float:
    data = graph[int(u)][int(v)]
    if not weight:
        return 1.0
    return float(data.get(weight, 1.0) or 1.0)


def connected_components(graph: Graph) -> Iterator[Set[int]]:
    visited: Set[int] = set()
    for start in graph.nodes():
        if start in visited:
            continue
        comp: Set[int] = set()
        queue = deque([int(start)])
        visited.add(int(start))
        while queue:
            node = int(queue.popleft())
            comp.add(node)
            for nbr in graph.neighbors(node):
                inbr = int(nbr)
                if inbr in visited:
                    continue
                visited.add(inbr)
                queue.append(inbr)
        yield comp


def number_connected_components(graph: Graph) -> int:
    return sum(1 for _ in connected_components(graph))


def is_connected(graph: Graph) -> bool:
    nodes = graph.nodes()
    if not nodes:
        return False
    return len(next(connected_components(graph))) == len(nodes)


def node_connected_component(graph: Graph, node: int) -> Set[int]:
    inode = int(node)
    if inode not in graph:
        return set()
    queue = deque([inode])
    comp = {inode}
    while queue:
        cur = int(queue.popleft())
        for nbr in graph.neighbors(cur):
            inbr = int(nbr)
            if inbr in comp:
                continue
            comp.add(inbr)
            queue.append(inbr)
    return comp


def _shortest_path_tree(
    graph: Graph,
    source: int,
    *,
    weight: Optional[str],
) -> Tuple[Dict[int, float], Dict[int, int]]:
    src = int(source)
    if src not in graph:
        return {}, {}
    if weight:
        dists: Dict[int, float] = {src: 0.0}
        prev: Dict[int, int] = {}
        heap: List[Tuple[float, int]] = [(0.0, src)]
        while heap:
            dist_u, u = heapq.heappop(heap)
            if dist_u > dists.get(u, float("inf")) + 1e-12:
                continue
            for v in graph.neighbors(u):
                iv = int(v)
                alt = float(dist_u + _edge_weight(graph, u, iv, weight))
                if alt + 1e-12 < dists.get(iv, float("inf")):
                    dists[iv] = alt
                    prev[iv] = u
                    heapq.heappush(heap, (alt, iv))
        return dists, prev

    dists = {src: 0.0}
    prev = {}
    queue = deque([src])
    while queue:
        u = int(queue.popleft())
        for v in graph.neighbors(u):
            iv = int(v)
            if iv in dists:
                continue
            dists[iv] = float(dists[u] + 1.0)
            prev[iv] = u
            queue.append(iv)
    return dists, prev


def has_path(graph: Graph, source: int, target: int) -> bool:
    dists, _ = _shortest_path_tree(graph, int(source), weight=None)
    return int(target) in dists


def shortest_path_length(
    graph: Graph,
    source: Optional[int] = None,
    target: Optional[int] = None,
    weight: Optional[str] = None,
):
    if source is None and target is None:
        def _iter_all_pairs():
            for src in graph.nodes():
                dists, _ = _shortest_path_tree(graph, int(src), weight=weight)
                yield int(src), {int(dst): float(val) for dst, val in dists.items()}

        return _iter_all_pairs()

    if source is None:
        raise ValueError("source must be provided when target is provided")

    dists, _ = _shortest_path_tree(graph, int(source), weight=weight)
    if target is None:
        return {int(dst): float(val) for dst, val in dists.items()}
    itarget = int(target)
    if itarget not in dists:
        raise NetworkXNoPath(f"No path between {source} and {target}")
    return float(dists[itarget])


def shortest_path(graph: Graph, source: int, target: int, weight: Optional[str] = None) -> List[int]:
    isource = int(source)
    itarget = int(target)
    dists, prev = _shortest_path_tree(graph, isource, weight=weight)
    if itarget not in dists:
        raise NetworkXNoPath(f"No path between {source} and {target}")
    path = [itarget]
    cur = itarget
    while cur != isource:
        cur = int(prev[cur])
        path.append(cur)
    path.reverse()
    return path


def all_pairs_dijkstra_path_length(graph: Graph, weight: str = "weight"):
    for src in graph.nodes():
        dists, _ = _shortest_path_tree(graph, int(src), weight=weight)
        yield int(src), {int(dst): float(val) for dst, val in dists.items()}


def bridges(graph: Graph) -> Iterator[Tuple[int, int]]:
    if isinstance(graph, MultiGraph):
        multiplicity = graph.multiplicity
    else:
        multiplicity = lambda u, v: 1

    visited: Set[int] = set()
    tin: Dict[int, int] = {}
    low: Dict[int, int] = {}
    timer = 0
    result: List[Tuple[int, int]] = []

    def dfs(node: int, parent: int) -> None:
        nonlocal timer
        visited.add(node)
        tin[node] = timer
        low[node] = timer
        timer += 1
        for nbr in graph.neighbors(node):
            v = int(nbr)
            if v == parent:
                continue
            if v in visited:
                low[node] = min(low[node], tin[v])
                continue
            dfs(v, node)
            low[node] = min(low[node], low[v])
            if low[v] > tin[node] and multiplicity(node, v) <= 1:
                result.append((node, v) if node <= v else (v, node))

    for node in graph.nodes():
        inode = int(node)
        if inode not in visited:
            dfs(inode, -1)

    seen: Set[Tuple[int, int]] = set()
    for edge in result:
        if edge in seen:
            continue
        seen.add(edge)
        yield edge


def k_edge_components(graph: Graph, k: int = 2) -> Iterator[Set[int]]:
    if int(k) != 2:
        yield from connected_components(graph)
        return
    bridge_set = {tuple(sorted(edge)) for edge in bridges(graph)}
    temp = Graph()
    temp.add_nodes_from(graph.nodes())
    for u, v in graph.edges():
        if tuple(sorted((u, v))) in bridge_set:
            continue
        temp.add_edge(u, v, **dict(graph[u][v]))
    yield from connected_components(temp)


def diameter(graph: Graph) -> int:
    if graph.number_of_nodes() <= 0:
        return 0
    max_dist = 0
    for src in graph.nodes():
        dists, _ = _shortest_path_tree(graph, int(src), weight=None)
        if len(dists) != graph.number_of_nodes():
            raise NetworkXNoPath("Graph is disconnected")
        local_max = max(int(val) for val in dists.values()) if dists else 0
        max_dist = max(max_dist, local_max)
    return int(max_dist)


def laplacian_matrix(graph: Graph, weight: str = "weight") -> _MatrixWrapper:
    nodes = sorted(int(node) for node in graph.nodes())
    idx = {node: i for i, node in enumerate(nodes)}
    size = len(nodes)
    arr = np.zeros((size, size), dtype=float)
    for u, v in graph.edges():
        w = max(float(graph[u][v].get(weight, 1.0) or 1.0), 1e-9)
        iu = idx[int(u)]
        iv = idx[int(v)]
        arr[iu, iu] += w
        arr[iv, iv] += w
        arr[iu, iv] -= w
        arr[iv, iu] -= w
    return _MatrixWrapper(arr)


def algebraic_connectivity(graph: Graph, weight: str = "weight", method: Optional[str] = None) -> float:
    if graph.number_of_nodes() < 2 or not is_connected(graph):
        return 0.0
    arr = laplacian_matrix(graph, weight=weight).toarray()
    evals = np.linalg.eigvalsh(arr)
    if len(evals) < 2:
        return 0.0
    return float(max(0.0, evals[1]))
