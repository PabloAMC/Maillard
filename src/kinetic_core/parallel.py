"""
src/kinetic_core/parallel.py -- INDEPENDENT POTS INTEGRATED IN PARALLEL, WITH THE CALLER'S STATE (2026-10-29).

WHY THIS EXISTS
===============
One evaluation of a fit objective integrates every pot (37 for B60's rows, about 70 for the joint fit) one after
another, and the pots are independent. The profile of 2026-10-29 (strategy M2): 3.9 s per evaluation on one core
with the inline pH coupling. Integrating the pots in parallel is the largest speed-up that cannot change a number.

THE ONE TRAP, AND HOW IT IS AVOIDED
===================================
The engine's switches are module globals (``sulfur.PH_COUPLING``, ``parameters_sulfur.CYS_AMADORI_ROUTE``,
``ph_state.B51``, a generator's ``install_reading`` ...), and some callers change them at every evaluation (the
route calibration does). A worker started by ``spawn`` imports the modules afresh and sees the DEFAULTS: it would
return wrong answers without any error. A long-lived pool would hold whatever state it was forked with. So the
pool is FORKED AT EACH CALL: every worker is a copy of the caller at that moment, switches included, and is gone
when the call returns. Forking 6 workers costs a few tenths of a second, small against a multi-second evaluation.
Results are bit-identical to the serial loop (tests/unit/test_parallel.py).

Workers: ``MAILLARD_WORKERS`` (default 1, i.e. serial, so nothing changes unless asked); never more than the
caller asks for. Where ``fork`` is unavailable the call runs serially, and so does a call made inside a worker.

WHERE IT IS USED, AND WHAT IT BUYS (measured 2026-10-29, min of three interleaved runs, 6 workers)
====================================================================================================
B2.3's ``run_systems`` / ``observables`` (and B64's profile variant): the pots of one evaluation, 2.5x with the
inline coupling (1.74 -> 0.69 s), 2.1x with the shipped fixed point (21.3 -> 10.3 s); the slowest pot (4 s) and
two pool forks per evaluation are the floor. ``fitting``'s Jacobians: the finite-difference points, which are
equal-sized whole evaluations, the coarse grain a fit spends most of its time on.
"""
from __future__ import annotations

import multiprocessing as mp
import os
from typing import Any, Callable, Iterable, List, Optional

_FN: Optional[Callable[[Any], Any]] = None
#: True inside a worker. A worker is a daemon and may not fork its own pool, so a nested call (a Jacobian's
#: points in parallel, each integrating its pots through ``map_ordered`` again) runs serially there.
_IN_WORKER = False


def workers(default: int = 1) -> int:
    """The worker count: ``MAILLARD_WORKERS`` if set, else ``default``."""
    try:
        return max(1, int(os.environ.get("MAILLARD_WORKERS", default)))
    except ValueError:
        return max(1, default)


def _mark_worker():
    global _IN_WORKER
    _IN_WORKER = True


def _call(item):
    return _FN(item)


def map_ordered(fn: Callable[[Any], Any], items: Iterable[Any], n_workers: Optional[int] = None) -> List[Any]:
    """``[fn(item) for item in items]``, computed by ``n_workers`` forked copies of this process when more than
    one is asked for and there is more than one item. Order is preserved; an exception in a worker is re-raised."""
    global _FN
    items = list(items)
    n = workers() if n_workers is None else max(1, int(n_workers))
    if n == 1 or len(items) < 2 or _IN_WORKER or "fork" not in mp.get_all_start_methods():
        return [fn(item) for item in items]
    _FN = fn  # read by the forked children through ``_call``; a closure needs no pickling this way
    try:
        with mp.get_context("fork").Pool(min(n, len(items)), initializer=_mark_worker) as pool:
            return pool.map(_call, items, chunksize=1)
    finally:
        _FN = None
