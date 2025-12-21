from __future__ import annotations

import heapq
from collections import defaultdict
from dataclasses import dataclass
from typing import Callable, DefaultDict, Dict, Generic, Iterable, List, Optional, Sequence, Tuple, TypeVar

from .linkscape_engine import UnionFind

T = TypeVar("T")


@dataclass(frozen=True)
class CircuitPick(Generic[T]):
    candidate: T
    corr_type: str
    score: float
    overlap_ratio: float


def select_circuit_utility(
    candidates: Iterable[T],
    *,
    budget: float,
    get_patch_ids: Callable[[T], Sequence[int]],
    get_pair_key: Callable[[T], Tuple[int, int]],
    get_cost: Callable[[T], float],
    get_base_roi: Callable[[T], float],
    overlap_ratio: Callable[[T, Sequence[object]], float],
    overlap_obj: Callable[[T], object],
    overlap_reject_ratio: float = 0.30,
    max_prior_per_pair: int = 3,
    diminishing_base: float = 0.5,
) -> Tuple[List[CircuitPick[T]], Dict[str, object]]:
    """
    Shared Circuit Theory greedy selector used by raster + vector.

    score = base_roi * multiplier

    multiplier:
      1.0  if candidate connects different components
      diminishing_base^(k+1) if already connected but spatially distinct (k = count already selected for this pair)
      0.01 if overlaps existing corridor for the pair by > overlap_reject_ratio
    """
    remaining = float(budget or 0.0)
    budget_used = 0.0
    primary_links = 0
    redundant_links = 0
    wasteful_links = 0

    uf = UnionFind()

    # Selected overlap representations by pair (limit memory per pair).
    selected_overlap_by_pair: DefaultDict[Tuple[int, int], List[object]] = defaultdict(list)
    selected_count_by_pair: DefaultDict[Tuple[int, int], int] = defaultdict(int)

    # Heap entries store a stamp so we can "re-score" candidates dynamically.
    heap: List[Tuple[float, int, int, T]] = []
    counter = 0
    stamps: Dict[int, int] = {}

    # Convert candidates to list once if iterable is one-shot.
    cand_list = list(candidates)

    def _stamp_key(cand: T) -> int:
        return id(cand)

    for cand in cand_list:
        base = float(get_base_roi(cand) or 0.0)
        cost = float(get_cost(cand) or 0.0)
        if base <= 0.0 or cost <= 0.0:
            continue
        counter += 1
        k = _stamp_key(cand)
        stamps[k] = stamps.get(k, 0) + 1
        heapq.heappush(heap, (-base, counter, stamps[k], cand))

    picks: List[CircuitPick[T]] = []

    while heap and remaining > 0:
        neg_score, _idx, stamp, cand = heapq.heappop(heap)
        k = _stamp_key(cand)
        if stamps.get(k, 0) != int(stamp):
            continue
        old_score = -float(neg_score)

        cost = float(get_cost(cand) or 0.0)
        if cost <= 0.0 or cost > remaining:
            continue

        base_roi = float(get_base_roi(cand) or 0.0)
        if base_roi <= 0.0:
            continue

        pids = [int(pid) for pid in get_patch_ids(cand) if pid is not None]
        if len(pids) < 2:
            continue

        try:
            roots = {int(uf.find(int(pid))) for pid in pids}
        except Exception:
            roots = {int(uf.find(int(pids[0]))), int(uf.find(int(pids[-1])))}

        pair = get_pair_key(cand)
        overlap_r = 0.0
        if len(roots) > 1:
            mult = 1.0
        else:
            prior = selected_overlap_by_pair.get(pair, [])
            overlap_r = float(overlap_ratio(cand, prior) or 0.0)
            if overlap_r > float(overlap_reject_ratio):
                mult = 0.01
            else:
                mult = float(diminishing_base) ** (int(selected_count_by_pair.get(pair, 0)) + 1)

        new_score = base_roi * mult
        if abs(new_score - old_score) > 1e-12:
            counter += 1
            stamps[k] = stamps.get(k, 0) + 1
            heapq.heappush(heap, (-new_score, counter, stamps[k], cand))
            continue

        if mult >= 1.0:
            corr_type = "primary"
            primary_links += 1
        else:
            corr_type = "redundant"
            redundant_links += 1
            if mult <= 0.01:
                corr_type = "wasteful"
                wasteful_links += 1

        # Commit: union all patch ids.
        for pid in pids:
            uf.find(int(pid))
        anchor = int(pids[0])
        for other in pids[1:]:
            uf.union(anchor, int(other))

        remaining -= cost
        budget_used += cost

        picks.append(CircuitPick(candidate=cand, corr_type=corr_type, score=float(new_score), overlap_ratio=float(overlap_r)))

        try:
            obj = overlap_obj(cand)
            if obj is not None:
                lst = selected_overlap_by_pair[pair]
                lst.append(obj)
                if len(lst) > int(max_prior_per_pair):
                    del lst[0]
        except Exception:
            pass
        selected_count_by_pair[pair] += 1

    stats: Dict[str, object] = {
        "budget_used": budget_used,
        "budget_remaining": remaining,
        "primary_links": primary_links,
        "redundant_links": redundant_links,
        "wasteful_links": wasteful_links,
        "corridors_used": len(picks),
    }
    return picks, stats
