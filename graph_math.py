"""
Graph Math Library for TerraLink
--------------------------------
Provides core graph-theoretic calculations for summary metrics such as entropy,
robustness, and redundancy.
"""

import math
import networkx as nx
import numpy as np
from typing import Dict, List, Tuple, Any

def build_graph_from_corridors(patches: Dict[int, Dict], corridors: List[Dict]) -> nx.Graph:
    """Helper to convert TerraLink patches/corridors to a NetworkX graph."""
    G = nx.Graph()
    for pid in patches:
        G.add_node(pid)
    
    for c in corridors:
        p1, p2 = c['patch1'], c['patch2']
        # Use distance if available, otherwise estimate from area (sqrt) or generic 1.0
        dist = c.get('distance_m', 1.0)
        G.add_edge(p1, p2, weight=dist, distance=dist, id=c.get('id'))
    return G

# --- ENTROPY CALCULATIONS ---

def calculate_movement_entropy(G: nx.Graph, alpha: float = 0.002) -> float:
    """
    Calculate movement entropy H_mov based on dispersal kernel.
    H_mov = -Σ p_i log(p_i)
    """
    H_total = 0.0
    n = G.number_of_nodes()
    if n == 0: return 0.0

    for node in G.nodes():
        neighbors = list(G.neighbors(node))
        if not neighbors:
            continue
            
        probs = []
        for nbr in neighbors:
            d = G[node][nbr].get('distance', 1.0)
            # Simplified dispersal kernel: exp(-alpha * distance)
            p = math.exp(-alpha * d)
            probs.append(p)
        
        # Normalize probabilities
        prob_sum = sum(probs)
        if prob_sum > 0:
            norm_probs = [p / prob_sum for p in probs]
            # Shannon entropy for this node
            H_i = -sum(p * math.log(p + 1e-10) for p in norm_probs)
            H_total += H_i

    return H_total / n

def calculate_topology_penalty(G: nx.Graph) -> float:
    """
    F(L): Penalizes deviation from dendritic (tree) structure.
    Returns number of cycles (extra edges beyond MST).
    """
    if G.number_of_nodes() == 0: return 0.0
    
    n_edges = G.number_of_edges()
    n_nodes = G.number_of_nodes()
    n_components = nx.number_connected_components(G)
    
    # For a forest (collection of trees), edges = nodes - components
    # Any excess edges imply cycles
    min_edges_needed = n_nodes - n_components
    cycles = max(0, n_edges - min_edges_needed)
    return float(cycles)

def calculate_disturbance_penalty(G: nx.Graph) -> float:
    """
    D(L): Penalizes long disturbance response times (diameter).
    Normalized by network size.
    """
    n = G.number_of_nodes()
    if n < 2: return 0.0
    
    if nx.is_connected(G):
        # Calculate diameter (longest shortest path)
        try:
            d = nx.diameter(G)
            return d / n
        except Exception:
            return 1.0
    else:
        # Penalize disconnection
        return 1.0

def calculate_total_entropy(G: nx.Graph, 
                          lambda_c: float = 1.0, 
                          lambda_f: float = 1.0, 
                          lambda_d: float = 1.0) -> Dict[str, float]:
    """Aggregate entropy metrics."""
    h_mov = calculate_movement_entropy(G)
    
    # Connectivity Constraint C(L): Simplified as number of components penalty
    n_comp = nx.number_connected_components(G)
    c_val = float(n_comp - 1) if n_comp > 0 else 0.0
    
    f_val = calculate_topology_penalty(G)
    d_val = calculate_disturbance_penalty(G)
    
    h_total = h_mov + (lambda_c * c_val) + (lambda_f * f_val) + (lambda_d * d_val)
    
    return {
        "H_total": h_total,
        "H_mov": h_mov,
        "C_connectivity": c_val,
        "F_topology": f_val,
        "D_disturbance": d_val
    }

# --- ROBUSTNESS CALCULATIONS ---

def calculate_two_edge_connectivity(G: nx.Graph) -> float:
    """
    Calculate rho_2: Fraction of node pairs with 2 edge-disjoint paths.
    Using a sampling approach for performance on large graphs, or full calculation for small ones.
    """
    n = G.number_of_nodes()
    if n < 3: return 0.0
    
    # For performance, we calculate the size of the 2-edge-connected components
    # This is a good proxy for rho_2 and much faster than all-pairs checks
    
    try:
        # k_edge_components returns sets of nodes in each 2-edge-connected component
        components = list(nx.k_edge_components(G, k=2))
        
        # Count pairs within these components
        two_connected_pairs = 0
        for comp in components:
            size = len(comp)
            two_connected_pairs += size * (size - 1) // 2
            
        total_pairs = n * (n - 1) // 2
        return two_connected_pairs / total_pairs if total_pairs > 0 else 0.0
    except Exception:
        return 0.0

def calculate_failure_probability(G: nx.Graph, k_failures: int = 1, iterations: int = 100) -> float:
    """
    Estimate P_fail(k): Probability that removing k edges disconnects the graph.
    Monte Carlo simulation.
    """
    edges = list(G.edges())
    m = len(edges)
    if m <= k_failures: return 1.0
    
    failures = 0
    for _ in range(iterations):
        # Randomly remove k edges
        G_test = G.copy()
        remove_indices = np.random.choice(m, k_failures, replace=False)
        for idx in remove_indices:
            G_test.remove_edge(*edges[idx])
            
        if not nx.is_connected(G_test) and nx.number_connected_components(G_test) > nx.number_connected_components(G):
             failures += 1
             
    return failures / iterations

def score_edge_for_loops(G_base: nx.Graph, u, v, weight) -> float:
    """
    Score a potential loop edge based on Betweenness Centrality Reduction.
    Higher score = better candidate for strategic redundancy.
    """
    # Simple heuristic: Connect nodes that are far apart in the current tree
    # This reduces diameter and provides shortcuts (high betweenness reduction)
    try:
        path_len = nx.shortest_path_length(G_base, u, v, weight='weight')
        # Score is path length in tree (shortcut value) divided by cost of new edge
        return path_len / (weight + 1e-6)
    except nx.NetworkXNoPath:
        return 0.0
