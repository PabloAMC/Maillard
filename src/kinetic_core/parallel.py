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

A JOB MAY CARRY ITS CONFIGURATION (A3 phase 1, 2026-11-01)
==========================================================
``map_ordered(fn, items, config=cfg)``, or an item wrapped as ``Configured(item, cfg)``: the worker installs the
``src.kinetic_core.config.Config`` at the start of the task (and restores its own state after it), so the answer no
longer depends on what the worker inherited. A spawn-started worker given the Config reproduces the forked one
(tests/unit/test_engine_config.py). Without a Config nothing changes: the pool is still forked at each call, and the
forked state is still what a job computes with. Script-level globals (a generator's own module state) are not in a
Config; a job that reads them still needs the fork.

Workers: ``MAILLARD_WORKERS`` (default 1, i.e. serial, so nothing changes unless asked); never more than the
caller asks for. Where ``fork`` is unavailable the call runs serially, and so does a call made inside a worker.

A PERSISTENT POOL, BEHIND A FLAG (A3 phase 2, 2026-11-01)
==========================================================
``MAILLARD_POOL=persistent`` (default ``fork_per_call``: nothing changes unless asked) keeps ONE pool of SPAWNED
workers for the life of the process, so a fit no longer pays a fork per evaluation. A spawned worker inherits
nothing, so every task travels as ``(fn, Configured(item, config))``: the item's own Config, else the caller's
configuration read at the call (``Config.from_globals()``), which the worker installs for the task and takes off
after it. Bit-identical to the forked pool (tests/unit/test_engine_config.py). Its limits:

* a function or item that cannot be pickled (a closure, a lambda: B2.3's and B64's per-system integrators, the
  signature scorers' draws) runs on the per-call fork, which carries it. Once the persistent pool exists the
  process has its handler threads, and Python warns that forking a multi-threaded process may deadlock; the child
  uses none of their locks, but that is why the flag stays off by default;
* state that is not an engine switch -- a script's own module globals -- is imported fresh in a spawned worker,
  so a job must carry it (the sub-model fits build each task's Config in the caller: ``configured``).

WHERE IT IS USED, AND WHAT IT BUYS (measured 2026-10-29, min of three interleaved runs, 6 workers)
====================================================================================================
B2.3's ``run_systems`` / ``observables`` (and B64's profile variant): the pots of one evaluation, 2.5x with the
inline coupling (1.74 -> 0.69 s), 2.1x with the shipped fixed point (21.3 -> 10.3 s); the slowest pot (4 s) and
two pool forks per evaluation are the floor. ``fitting``'s Jacobians: the finite-difference points, which are
equal-sized whole evaluations, the coarse grain a fit spends most of its time on.
"""
from __future__ import annotations

import atexit
import multiprocessing as mp
import os
import pickle
from dataclasses import dataclass
from typing import Any, Callable, Iterable, List, Optional

_FN: Optional[Callable[[Any], Any]] = None
#: True inside a worker. A worker is a daemon and may not fork its own pool, so a nested call (a Jacobian's
#: points in parallel, each integrating its pots through ``map_ordered`` again) runs serially there.
_IN_WORKER = False
#: the persistent pool, ``(size, pool)``, made on first use when ``MAILLARD_POOL=persistent``
_PERSISTENT = None
POOL_MODES = ("fork_per_call", "persistent")


def workers(default: int = 1) -> int:
    """The worker count: ``MAILLARD_WORKERS`` if set, else ``default``."""
    try:
        return max(1, int(os.environ.get("MAILLARD_WORKERS", default)))
    except ValueError:
        return max(1, default)


def pool_mode() -> str:
    """``MAILLARD_POOL``: ``fork_per_call`` (the default) or ``persistent``; anything else is an error."""
    mode = os.environ.get("MAILLARD_POOL", "fork_per_call")
    if mode not in POOL_MODES:
        raise ValueError(f"MAILLARD_POOL={mode!r}: one of {POOL_MODES}")
    return mode


def _mark_worker():
    global _IN_WORKER
    _IN_WORKER = True


@dataclass(frozen=True)
class Configured:
    """An item and the engine ``Config`` it is computed under (``src.kinetic_core.config``)."""

    item: Any
    config: Any


def run_job(fn: Callable[[Any], Any], job: Any) -> Any:
    """``fn(job)``; for a ``Configured`` job, ``fn(job.item)`` with its Config installed for the task."""
    if isinstance(job, Configured):
        from .config import applied

        with applied(job.config):
            return fn(job.item)
    return fn(job)


def run_pickled(payload) -> Any:
    """``run_job`` for a pool that pickles its tasks (spawn): ``payload`` is ``(fn, job)``."""
    fn, job = payload
    return run_job(fn, job)


def configured(items: Iterable[Any], config_of: Optional[Callable[[Any, Any], Any]] = None) -> List[Configured]:
    """Every item with the Config it runs in, built HERE from this process's configuration at this moment (what a
    worker forked now would inherit): ``config_of(item, base)``, or ``base`` itself. A job built this way computes
    the same answer in any worker, forked, spawned or persistent."""
    from .config import Config

    base = Config.from_globals()
    return [Configured(item, base if config_of is None else config_of(item, base)) for item in items]


def run_under(job: Configured, fn: Callable[[], Any]) -> Any:
    """``fn()`` in this process with ``job``'s Config installed, then the state back: to read a job's result in the
    configuration it was computed in (an observation operator reads the pH switches)."""
    from .config import applied

    with applied(job.config):
        return fn()


def _call(item):
    return run_job(_FN, item)


def _persistent_pool(n: int):
    global _PERSISTENT
    if _PERSISTENT is not None and _PERSISTENT[0] != n:
        shutdown_persistent_pool()
    if _PERSISTENT is None:
        _PERSISTENT = (n, mp.get_context("spawn").Pool(n, initializer=_mark_worker))
    return _PERSISTENT[1]


def shutdown_persistent_pool() -> None:
    """Close the persistent pool, if there is one (also at exit)."""
    global _PERSISTENT
    if _PERSISTENT is not None:
        _, pool = _PERSISTENT
        _PERSISTENT = None
        pool.terminate()
        pool.join()


atexit.register(shutdown_persistent_pool)


def _persistent_payloads(fn, items) -> Optional[list]:
    """``[(fn, Configured job)]`` for the persistent pool, every job carrying a Config (the caller's where the item
    has none), or None when the function or an item cannot be pickled (then the per-call fork runs it)."""
    caller = None
    jobs = []
    for item in items:
        if not isinstance(item, Configured):
            if caller is None:
                from .config import Config

                caller = Config.from_globals()
            item = Configured(item, caller)
        jobs.append(item)
    payloads = [(fn, job) for job in jobs]
    try:
        pickle.dumps(payloads)
    except (pickle.PicklingError, TypeError, AttributeError):
        return None
    return payloads


def map_ordered(fn: Callable[[Any], Any], items: Iterable[Any], n_workers: Optional[int] = None, *,
                config: Any = None) -> List[Any]:
    """``[fn(item) for item in items]``, computed by ``n_workers`` forked copies of this process when more than
    one is asked for and there is more than one item. Order is preserved; an exception in a worker is re-raised.
    ``config`` (a ``Config``) is installed for every item that does not carry its own (``Configured``). With
    ``MAILLARD_POOL=persistent`` the workers are a persistent spawned pool and every job carries a Config."""
    global _FN
    items = list(items)
    if config is not None:
        items = [item if isinstance(item, Configured) else Configured(item, config) for item in items]
    n = workers() if n_workers is None else max(1, int(n_workers))
    if n == 1 or len(items) < 2 or _IN_WORKER or "fork" not in mp.get_all_start_methods():
        return [run_job(fn, item) for item in items]
    if pool_mode() == "persistent":
        payloads = _persistent_payloads(fn, items)
        if payloads is not None:
            return _persistent_pool(n).map(run_pickled, payloads, chunksize=1)
    _FN = fn  # read by the forked children through ``_call``; a closure needs no pickling this way
    try:
        with mp.get_context("fork").Pool(min(n, len(items)), initializer=_mark_worker) as pool:
            return pool.map(_call, items, chunksize=1)
    finally:
        _FN = None
