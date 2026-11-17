# SPDX-License-Identifier: LGPL-3.0-or-later
"""``axtb.program``: axtb as a program of qccompute's protocol, structurally (no qccompute
import): its name, version, capabilities, schemas, plans and calculations."""

from __future__ import annotations

from typing import Any

from qcdata import CalcType, ProgramInput, ProgramOutput, SinglePointData

from ._keywords import CALCTYPES, METHODS, plan_input, plan_threads, settings_schema
from ._qc import VERSION, compute

__all__ = ["Axtb", "program"]


_RECORD_SCHEMA: dict[str, Any] = {
    "$schema": "https://json-schema.org/draft/2020-12/schema",
    "title": "axtb record",
    "type": "object",
    "description": "axtb's detail of a calculation; diagnostic, not stable across versions.",
    "properties": {
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
                "attempts": {"type": "array", "items": {"type": "object"}},
            },
        },
        "energy_subtotals": {
            "type": "object",
            "description": "Terms that overlap the energy_components (GFN2-xTB's and g-xTB's electronic energy).",
            "additionalProperties": {"type": "number"},
        },
        "contingent": {
            "type": "object",
            "description": "Contingent settings and whether each fired (key -> {fired}).",
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
            "sessions": False,
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


#: The object registered with qccompute.
program = Axtb()
