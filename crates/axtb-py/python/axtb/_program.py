# SPDX-License-Identifier: LGPL-3.0-or-later
"""``axtb.program``: axtb as a program of qccompute's protocol, structurally (no qccompute
import): its name, version, capabilities, schemas, plans and calculations."""

from __future__ import annotations

from typing import Any

from qcdata import CalcType, ProgramInput, ProgramOutput, SinglePointData

from . import _axtb
from ._keywords import CALCTYPES, METHODS, alias, canonical_method, plan_input, plan_threads, settings_schema
from ._qc import VERSION, compute
from ._session import Session

__all__ = ["Axtb", "program"]


_RECORD_SCHEMA: dict[str, Any] = {
    "$schema": "https://json-schema.org/draft/2020-12/schema",
    "title": "axtb record",
    "type": "object",
    "description": "axtb's detail of a calculation; diagnostic, not stable across versions.",
    "properties": {
        "warm_start": {
            "type": ["boolean", "null"],
            "description": "Session frames: whether the SCF started from its lineage's history (null: GFN-FF, which "
            "keeps none). The composite: whether history was installed and a part's restart policy reuses it.",
        },
        "from_lineage": {"description": "Session frames: the lineage (or seed) whose history started it, or null."},
        "accuracy": {
            "type": "object",
            "description": "With an accuracy target: the target, the thresholds used and the model's error "
            "estimate at them (not measured on this run).",
            "properties": {
                "target": {"type": "object"},
                "thresholds": {"type": "object"},
                "estimate": {"type": "object"},
            },
        },
        "numerics": {"type": "string", "description": "The numerics identity of the code that made the numbers."},
        "method": {"type": "string", "description": "The canonical method."},
        "conventions": {"type": "string", "description": "The canonical conventions."},
        "eigensolver": {"type": "string", "description": "faer, or the run-time LAPACK library (AXTB_LAPACK)."},
        "scf": {
            "type": "object",
            "description": "How the SCF went.",
            "properties": {
                "start": {"enum": ["eeq", "sad", "eeq-bc", "previous", "given"]},
                "iterations": {"type": "integer"},
                "mixer_used": {"type": "string"},
                "attempts": {
                    "type": "array",
                    "items": {"type": "object"},
                    "description": "GFN2-xTB's mixer attempts; with scf.retries every try of the ladder (retry, "
                    "etemp, max_iter, converged, iterations, start or error; the composite's part), in a failed "
                    "output too.",
                },
                "restart_fallback": {
                    "type": "boolean",
                    "description": "Session frames: the SCF from history failed and ran again fresh.",
                },
            },
        },
        "energy_subtotals": {
            "type": "object",
            "description": "Terms that overlap the energy_components (GFN2-xTB's and g-xTB's electronic energy).",
            "additionalProperties": {"type": "number"},
        },
        "contingent": {
            "type": "object",
            "description": "Contingent settings and whether each fired (key -> {fired}; scf.retries: and the "
            "rung that gave the result, the composite's parts that climbed).",
            "additionalProperties": {"type": "object", "properties": {"fired": {"type": "boolean"}}},
        },
        "hessian": {
            "type": "object",
            "properties": {"method": {"type": "string"}, "max_asymmetry": {"type": "number"}},
        },
        "composite": {
            "type": "object",
            "description": "The g-xTB solvation composite: definition, energies of its parts, their SCFs.",
        },
        "wavefunction": {
            "type": "object",
            "description": "Which part's wavefunction the outputs are (the composite: gxtb).",
        },
        "gfnff": {
            "type": "object",
            "description": "GFN-FF: hydrogen and halogen bonds, the fragment charges and the topology.",
        },
        "error": {
            "type": "object",
            "description": "A failure's fields: z of an unsupported element; iterations, energy_change and residual "
            "of an SCF that did not converge.",
        },
    },
}


class Axtb:
    """axtb as a qccompute program (``axtb.program``): the object registered under the
    ``qccompute.programs`` entry points. It needs no qccompute: capabilities and plans are
    plain dicts, which qccompute validates."""

    name = "axtb"

    def __repr__(self) -> str:
        return "axtb.program"

    def version(self) -> str:
        """axtb's version (PEP 440)."""
        return VERSION

    def capabilities(self) -> dict[str, Any]:
        """Energies, gradients and Hessians of every method; in process and thread-safe."""
        calctypes = sorted(CALCTYPES)
        return {
            "calctypes": calctypes,
            "methods": {m: calctypes for m in METHODS},
            "sessions": True,
            "batches": False,
            "in_process": True,
            "thread_safe": True,
        }

    def settings_schema(self, calctype: CalcType | str) -> dict[str, Any]:
        """The JSON Schema of the keywords of a calctype."""
        return settings_schema(CalcType(calctype).value)

    def record_schema(self, calctype: CalcType | str) -> dict[str, Any]:
        """The JSON Schema of ``record``."""
        return _RECORD_SCHEMA

    def plan(self, inp: ProgramInput) -> dict[str, Any]:
        """Check an input and resolve its settings, without running it: ``resolved`` (threads
        deferred to the run unless the keyword gives them), ``resources`` (the threads the
        calculation asks for: see ``docs/python.md``), ``warnings`` and ``record``.

        Raises:
            ProgramError: The input is refused (kind ``input`` or ``unsupported``).
        """
        p = plan_input(inp)
        settings = dict(p.settings)
        deferred: list[str] = []
        if p.threads is not None:
            settings["threads"] = p.threads
        else:
            deferred.append("threads")
        return {
            "resolved": {"settings": settings, "reasons": p.reasons, "deferred": deferred},
            "resources": {"threads": plan_threads(p, inp.calctype.value)},
            "warnings": [],
            "record": p.record,
        }

    def compute(self, inp: ProgramInput, ctx: Any = None) -> ProgramOutput[ProgramInput, SinglePointData]:
        """`axtb.compute`."""
        return compute(inp, ctx)

    def open_session(self, inp: ProgramInput, ctx: Any = None, calctypes: Any = None) -> Session:
        """`axtb.Session`: the protocol's session object."""
        return Session(inp, calctypes=calctypes, ctx=ctx)

    def numerics(self, spec: Any) -> str:
        """The numerics identity of a spec (an input or a subprogram's `ProgramSpec`): the method's
        epoch (``gfn2:1``; the g-xTB solvation composite ``gxtb-v2.0.1:1+gfn2:1``), the same for every
        conventions of a method; axtb's version for a spec that names no method axtb has (its
        result is a refusal). Never raises: the runner calls it before looking up the store."""
        try:
            method = canonical_method(spec.model.method)
            solvent = (spec.keywords or {}).get("solvent") or {}
            composite = (
                method.startswith("gxtb")
                and isinstance(solvent, dict)
                and isinstance(solvent.get("correction"), str)
                and alias(solvent["correction"]) in ("gfn2", "gfn2xtb")
                and solvent.get("name") is not None
            )
            return str(_axtb.numerics(method, composite))
        except Exception:  # noqa: BLE001 - an input axtb refuses: its result follows axtb's checks
            return f"axtb {self.version()}"


#: The object registered with qccompute.
program = Axtb()
