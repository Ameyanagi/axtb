# SPDX-License-Identifier: LGPL-3.0-or-later
"""``axtb.Session``: one molecule at many geometries on one engine, warm-started per lineage
(qccompute's session protocol v4; qcdata's `Request` and `Frame`).

The input pins the atoms, charge, multiplicity, method and keywords; its geometry and calctype
are a request's defaults. The engine (`axtb._axtb.Engine`) builds the method once (GFN-FF its
topology too) and keeps one committed wavefunction per lineage: a request starts from its
lineage's history, else from the first of its seeds with one, else fresh, as the method's restart
policy says; a success with ``commit`` becomes the lineage's history. A request runs on the
threads of its calctype and size, on the engine's calculator or a fork of it on those threads;
``compute_many`` runs different lineages side by side as their threads fit into the session's.
"""

from __future__ import annotations

import json
import threading
import time
from collections.abc import Iterable, Sequence
from typing import Any, NamedTuple

import numpy as np
from qcdata import (
    Accuracy,
    ExecutionInfo,
    Frame,
    Lineage,
    ProgramError,
    ProgramErrors,
    ProgramInput,
    ProgramWarning,
    Request,
    SinglePointData,
    Structure,
    lineage_key,
)

from . import _axtb, _qc
from ._errors import from_driver, refuse
from ._keywords import CALCTYPES, largest, plan_input
from ._models import build

__all__ = ["Session"]


def _key(lineage: Any) -> str:
    """`qcdata.lineage_key`, with the common keys fast."""
    if lineage is None:
        return "null"
    if type(lineage) is str and lineage.isalnum() and lineage.isascii():
        return f'"{lineage}"'
    return lineage_key(lineage)


_LOWER = {"energy": ("energy",), "gradient": ("energy", "gradient"), "hessian": ("energy", "gradient", "hessian")}
#: The outputs a request may name.
OUTPUTS = ("orbitals", "overlap", "hamiltonian", "core_hamiltonian", "bond_orders", "ir", "dipole_moment")


def session_calctypes(calctypes: Iterable[Any] | None) -> frozenset[str]:
    """The calctypes a session may run: every calctype (None), or each declared one with its lower
    derivatives (qccompute's ``session_calctypes``, without importing qccompute)."""
    if not calctypes:
        return frozenset(CALCTYPES)
    allowed: set[str] = set()
    for c in calctypes:
        name = getattr(c, "value", c)
        if name not in _LOWER:
            raise refuse(
                "unsupported",
                "unsupported_calctype",
                f"calctypes: axtb sessions run energy, gradient and hessian, not {name}",
            )
        allowed.update(_LOWER[name])
    return frozenset(allowed)


class _Prepared(NamedTuple):
    """A checked request: what the engine runs (``args``: positions, task, lineage key, seed
    keys, commit, SCF accuracy, orbitals, threads) and what its frame needs."""

    args: tuple[Any, ...]
    x: np.ndarray
    calctype: str
    wanted: frozenset[str]
    lineage: Any
    key: str
    seeds: list[Any]
    acc_record: dict[str, Any] | None
    acc_warnings: list[Any]


class Session:
    """A session of axtb: the protocol's ``Session`` (``axtb.program.open_session``) and axtb's
    own way to run one molecule at many geometries.

    Every calculation runs on its own threads: the ``threads`` keyword, else the plan's for its
    calctype and size, at most the grant. ``compute_many`` runs different lineages side by side as
    long as their threads fit into the session's threads in all: the grant, and without one (a
    direct call) this process's CPUs.

    Args:
        inp: The atoms, charge, multiplicity, method and keywords, pinned for the session's
            lifetime; its geometry and calctype are a request's defaults.
        calctypes: The calctypes it may run (each with its lower derivatives); None: energy,
            gradient and hessian. ``hessian.*`` and ``freq.*`` keywords need ``hessian`` among them.
        ctx: A runner's grant (`qccompute.program.RunContext`); only ``ctx.threads`` is read.

    Raises:
        ProgramError: The input is refused (kind ``input`` or ``unsupported``), or the method
            cannot be set up for the molecule (a GFN-FF topology).
        TypeError: ``inp`` is not a `qcdata.ProgramInput`.
    """

    def __init__(self, inp: ProgramInput, *, calctypes: Iterable[Any] | None = None, ctx: Any = None):
        if not isinstance(inp, ProgramInput):
            raise TypeError(f"axtb.Session opens from a qcdata.ProgramInput, not {type(inp).__name__}")
        allowed = session_calctypes(calctypes)
        p = plan_input(inp, calctypes=allowed, session=True)
        per, total, warnings = _qc.session_threads(p, allowed, ctx)
        # small molecules run on one thread (the same bits), as compute does
        used = {
            c: 1 if len(p.numbers) <= _axtb.single_thread_atoms(p.method, c == "hessian") else n
            for c, (n, _) in per.items()
        }
        # the engine's own calculator runs on the input calctype's threads
        own = inp.calctype.value if inp.calctype.value in allowed else largest(allowed)
        request = json.dumps({**p.request, "threads": used[own]}, separators=(",", ":"))
        try:
            self._engine = _axtb.Engine(request, p.numbers, p.positions)
            if p.start is not None:
                # every request without history starts from the given results
                self._engine.set_given(**_qc.given_args(p))
        except _axtb.Error as x:
            error, _ = from_driver(x)
            raise ProgramError(error) from None
        self.calctypes = allowed
        self.program = "axtb"
        self.program_version = _qc.VERSION
        self._inp = inp
        self._p = p
        self._structure = inp.require_structure()
        self._natoms = len(p.numbers)
        # each calctype's threads (and why), the threads it runs on, the input's calctype and the
        # session's threads in all
        self._threads = per
        self._used = used
        self._own = own
        self._total = total
        self._warnings = warnings
        self._lock = threading.Lock()
        self._tasks: dict[tuple[str, frozenset[str], float | None], str] = {}
        self._static: dict[str, Any] | None = None
        self._geometry = np.asarray(self._structure.geometry, dtype=float).reshape(-1, 3)
        self._default_outputs = _qc.implied(p.outputs)
        self._scale = float(p.settings.get("accuracy_scale", 1.0))
        self._closed = False

    # ------------------------------------------------------------------ the protocol

    @property
    def input_data(self) -> ProgramInput:
        """The input the session was opened with."""
        return self._inp

    def compute(
        self,
        positions: Any = None,
        calctype: Any = None,
        *,
        accuracy: Accuracy | float | None = None,
        lineage: Lineage = None,
        outputs: Iterable[str] = (),
        seed: Iterable[Lineage] = (),
        commit: bool = True,
    ) -> Frame:
        """Evaluate one geometry.

        Args:
            positions: Positions (bohr), n_atoms x 3; None: the input's geometry.
            calctype: ``energy``, ``gradient`` or ``hessian``; None: the input's.
            accuracy: Accuracy targets of this request (a bare number targets the calctype's
                highest derivative); None: the input's ``accuracy`` keyword, else
                ``accuracy_scale``.
            lineage: The lineage this request continues (any JSON value; None: the default).
            outputs: Optional outputs (``orbitals``, ``overlap``, ``hamiltonian``,
                ``bond_orders``; ``ir`` with a Hessian), beside the keywords' ``output.*``.
            seed: Lineages to start from, in order, when ``lineage`` has no history.
            commit: Whether a success becomes ``lineage``'s history.

        Returns:
            The frame: the structure evaluated, its results, warnings and record (``warm_start``,
            ``from_lineage``).

        Raises:
            ProgramError: The request is refused or failed (``exc.frame``: the failed frame); the
                session stays usable (a panic closes it) and the lineage's history is untouched.
        """
        start = time.perf_counter()
        with self._lock:
            q = self._prepare(positions, calctype, accuracy, lineage, outputs, seed, commit)
            r: dict[str, Any] | _axtb.Error
            try:
                r = self._engine.compute(*q.args)
            except _axtb.Error as e:
                r = e
            return self._finish(q, r, start)

    def compute_many(self, requests: Sequence[Request]) -> list[Frame]:
        """Evaluate many requests: each lineage's in order, each reading its lineage's history as
        it stands at its place, and different lineages side by side as long as their threads fit
        into the session's threads in all (`Session`), taken in request order. A request whose
        lineage has no history waits for its seeds' earlier requests. The frames are what
        one-at-a-time ``compute`` calls in request order give, bit for bit; failures stay in
        place.

        Returns:
            One frame per request, in order, all successful.

        Raises:
            ProgramErrors: Some requests failed; ``results`` holds every frame.
        """
        frames: list[Frame | None] = [None] * len(requests)
        with self._lock:
            ready: list[tuple[int, float, _Prepared]] = []
            for i, req in enumerate(requests):
                start = time.perf_counter()
                try:
                    q = self._prepare(
                        req.positions, req.calctype, req.accuracy, req.lineage, req.outputs, req.seed, req.commit
                    )
                except ProgramError as e:
                    frames[i] = e.frame if e.frame is not None else self._failed(req.positions, e, None)
                    continue
                ready.append((i, time.perf_counter() - start, q))
            if ready:
                out = self._engine.compute_many([q.args for *_, q in ready], self._total)
                for (i, checks, q), r in zip(ready, out, strict=True):
                    # the request's own time: its checks, its calculation and its frame
                    took = checks + (r["wall"] if isinstance(r, dict) else 0.0)
                    try:
                        frames[i] = self._finish(q, r, time.perf_counter() - took)
                    except ProgramError as e:
                        frames[i] = e.frame
        return ProgramErrors.check(frames)

    def reset(self, lineage: Lineage = ...) -> None:
        """Forget the history of every lineage (no argument), of the default lineage (None) or of
        ``lineage``. The topology and the settings stay."""
        with self._lock:
            if not self._closed:
                self._engine.reset(None if lineage is ... else lineage_key(lineage))

    def plan(self) -> dict[str, Any]:
        """The settings pinned at opening (with the threads of the input's calctype), the
        resources (the most threads a calculation takes), the warnings and the record (method,
        conventions, calctypes, numerics; ``threads_used`` by the input's calctype,
        ``threads_by_calctype``, and ``lineages_at_once``: how many requests of the input's
        calctype ``compute_many`` runs side by side)."""
        n, why = self._threads[self._own]
        used = self._used[self._own]
        return {
            "resolved": {
                "settings": {**self._p.settings, "threads": n},
                "reasons": {**self._p.reasons, "threads": why},
                "deferred": [],
            },
            "resources": {"threads": max(m for m, _ in self._threads.values())},
            "warnings": [w.model_dump() for w in self._warnings],
            "record": {
                **self._p.record,
                "threads_used": used,
                "threads_by_calctype": {c: self._used[c] for c in CALCTYPES if c in self._used},
                "lineages_at_once": max(1, self._total // used),
            },
        }

    def close(self) -> None:
        """Drop the engine and every lineage's history; later requests fail (``session_closed``)."""
        with self._lock:
            if not self._closed:
                self._engine.close()
                self._closed = True

    def lineages(self) -> list[Lineage]:
        """The lineages with history (their JSON values), beyond the protocol."""
        with self._lock:
            return [] if self._closed else [json.loads(k) for k in self._engine.lineages()]

    def __enter__(self) -> Session:
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()

    def __repr__(self) -> str:
        state = "closed" if self._closed else f"{len(self._engine.lineages())} lineages"
        return f"axtb.Session({self._p.method}, {self._natoms} atoms, {state})"

    # ------------------------------------------------------------------ one request

    def _structure_at(self, x: np.ndarray) -> Any:
        """The session's structure at positions `x` (its other fields shared, as model_copy does)."""
        return build(Structure, **{**self._structure.__dict__, "geometry": x})

    def _failed(
        self,
        positions: Any,
        e: ProgramError,
        x: np.ndarray | None,
        record: dict | None = None,
        threads: int | None = None,
    ) -> Frame:
        """The failed frame of a request: the structure evaluated, empty results, the error, the
        threads it ran on (default: the input calctype's)."""
        if x is None:
            x = np.asarray(self._structure.geometry if positions is None else positions, dtype=float)
            if x.size != 3 * self._natoms:
                x = np.asarray(self._structure.geometry, dtype=float)
            x = x.reshape(-1, 3)
        return build(
            Frame,
            structure=self._structure_at(x),
            results=build(SinglePointData, provenance=_qc.provenance()),
            success=False,
            traceback=e.error.message,
            error=e.error,
            record=record or {},
            execution=build(ExecutionInfo, threads=self._used[self._own] if threads is None else threads),
        )

    def _positions(self, positions: Any) -> np.ndarray:
        """The positions of a request (the engine checks their shape and that they are finite)."""
        if positions is None:
            return self._geometry
        x = np.asarray(positions, dtype=float)
        return x.reshape(self._natoms, 3) if x.size == 3 * self._natoms else x

    def _outputs(self, outputs: Iterable[str], calctype: str) -> frozenset[str]:
        p = self._p
        if not outputs:
            return self._default_outputs
        wanted = set(outputs)
        unknown = wanted - set(OUTPUTS)
        if unknown:
            raise refuse(
                "unsupported",
                "not_applicable",
                f"outputs: axtb gives no {', '.join(sorted(unknown))}; it gives {', '.join(OUTPUTS)}",
            )
        wanted.discard("dipole_moment")  # always there
        if "ir" in wanted:
            if calctype != "hessian":
                raise refuse("unsupported", "not_applicable", "outputs: ir comes with a hessian request")
            wanted.discard("ir")  # with the Hessian whenever its method gives dipole derivatives
        if "core_hamiltonian" in wanted:
            raise refuse("unsupported", "not_applicable", "outputs: core_hamiltonian: tight binding's H0 is model data")
        if wanted and not p.scf:
            raise refuse(
                "unsupported",
                "not_applicable",
                f"outputs: GFN-FF has no orbitals, SCF matrices or bond orders ({', '.join(sorted(wanted))})",
            )
        if "hamiltonian" in wanted and p.method == "gxtb-v2.0.1" and p.program_conventions:
            raise refuse(
                "unsupported",
                "not_applicable",
                "outputs: hamiltonian: g-xTB 2.0.1 in its program's conventions diagonalizes a mixed matrix last",
            )
        return _qc.implied(wanted | p.outputs)

    def _prepare(
        self,
        positions: Any,
        calctype: Any,
        accuracy: Any,
        lineage: Any,
        outputs: Iterable[str],
        seed: Iterable[Any],
        commit: bool,
    ) -> _Prepared:
        """Check a request and resolve what the engine runs; a refusal raises `ProgramError`
        with its failed frame."""
        p = self._p
        x = None
        try:
            if self._closed:
                raise refuse("input", "session_closed", "the session is closed")
            calctype = self._inp.calctype.value if calctype is None else getattr(calctype, "value", calctype)
            x = self._positions(positions)
            if calctype not in self.calctypes:
                raise refuse(
                    "input",
                    "not_applicable",
                    f"calctype: this session was opened for {sorted(self.calctypes, key=CALCTYPES.index)}, not "
                    f"{calctype}; open it with calctypes that include it",
                )
            wanted = self._outputs(outputs, calctype)
            target = Accuracy.of(accuracy, calctype) if accuracy is not None else p.accuracy
            acc_warnings: list[ProgramWarning]
            if target is None:
                factor, h, acc_record, acc_warnings = self._scale, None, None, []
            else:
                factor, h, acc_record, acc_warnings = _qc.accuracy(p, target, calctype)
            seeds = list(seed)
            key = _key(lineage)
            task = self._tasks.get((calctype, wanted, h))
            if task is None:
                task = self._tasks[(calctype, wanted, h)] = _qc.task(p, calctype, wanted, h)
        except ProgramError as e:
            if e.frame is None:
                e.frame = self._failed(positions, e, x if x is not None and x.shape == (self._natoms, 3) else None)
            raise
        orbitals = _qc.needs_orbitals(wanted)
        threads = self._used[calctype]
        args = (x, task, key, [_key(s) for s in seeds], commit, factor if p.scf else None, orbitals, threads)
        return _Prepared(args, x, calctype, wanted, lineage, key, seeds, acc_record, acc_warnings)

    def _finish(self, q: _Prepared, r: dict[str, Any] | _axtb.Error, start: float) -> Frame:
        """The frame of a request from its engine response; an engine error raises its
        `ProgramError`, with the failed frame."""
        if isinstance(r, _axtb.Error):
            error, fields = from_driver(r)
            if r.kind in ("panic", "closed"):
                self._closed = True
            exc = ProgramError(error)
            rec: dict[str, Any] = {"error": fields}
            if r.attempts:
                rec["scf"] = {"attempts": r.attempts}
            exc.frame = self._failed(None, exc, q.x, rec, q.args[-1])
            raise exc from None
        return self._frame(start, r, q)

    def _frame(self, start: float, r: dict[str, Any], q: _Prepared) -> Frame:
        """The frame of a request's engine response `r`; `start`: when the request began."""
        p = self._p
        x, calctype, wanted, lineage, key, seeds = q.x, q.calctype, q.wanted, q.lineage, q.key, q.seeds
        acc_record, acc_warnings = q.acc_record, q.acc_warnings
        if self._static is None:
            self._static = _qc.static_results(r, p)
        res = _qc.results(r, p, calctype, wanted, static=self._static)
        history = r["history"]
        from_lineage = None
        if history is not None:
            if history == key:
                from_lineage = lineage
            else:
                from_lineage = next(s for s in seeds if _key(s) == history)
        scf = r.get("scf")
        if scf is None:
            warm: bool | None = None
        elif p.composite:
            warm = history is not None and any(p.restart)
        else:
            warm = scf["guess"] == "previous"
        if r["given"]:
            warm = False
        rec: dict[str, Any] = {"warm_start": warm, "from_lineage": from_lineage}
        rec.update(_qc.record(r, p, session=True))
        if r["given"] and "scf" in rec:
            rec["scf"]["start"] = "given"
        if scf is not None and not p.composite and history is not None and p.restart[0] and not warm:
            rec["scf"]["restart_fallback"] = True
        contingent = rec.setdefault("contingent", {})
        if "scf.restart" in p.contingent:
            contingent["scf.restart"] = {"fired": bool(warm)}
        for path in p.contingent:
            if path.startswith(("hessian.", "freq.")):
                contingent[path] = {"fired": calctype == "hessian"}
        if not contingent:
            del rec["contingent"]
        if acc_record is not None:
            rec["accuracy"] = acc_record
        return build(
            Frame,
            structure=self._structure_at(x),
            results=res,
            success=True,
            warnings=acc_warnings + _qc.warnings_of(r, res),
            record=rec,
            execution=build(ExecutionInfo, threads=r["threads"], timings=_qc.timings(r, time.perf_counter() - start)),
        )
