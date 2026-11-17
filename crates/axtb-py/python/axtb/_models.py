# SPDX-License-Identifier: LGPL-3.0-or-later
"""Building qcdata models from trusted values, fast.

axtb's outputs are built from the engine's numbers, of the right types and shapes. Validating
them again costs more than a small molecule's whole calculation (a water gradient: ~0.2 ms in
the engine), so `build` makes models the way pydantic's ``model_construct`` does, without
validation, and with each class's defaults computed once (``model_construct`` computes every
default, and inspects every default factory's signature, on every call). The tests validate
what it builds (``model_validate`` of the dump equals the model, and dumps the same JSON).
"""

from __future__ import annotations

import copy
from typing import Any

from pydantic import BaseModel

__all__ = ["build"]

_REQUIRED, _FIXED, _EMPTY, _COPY = range(4)
# per class: the fields in order, each (name, kind, default), and whether model_post_init runs
_TEMPLATES: dict[type, tuple[list[tuple[str, int, Any]], bool]] = {}


def _template(cls: Any) -> tuple[list[tuple[str, int, Any]], bool]:
    fields = []
    for name, f in cls.model_fields.items():
        if f.is_required():
            fields.append((name, _REQUIRED, None))
            continue
        value = f.get_default(call_default_factory=True, validated_data={})
        if isinstance(value, dict | list | set) and not value:
            fields.append((name, _EMPTY, type(value)))
        elif isinstance(value, dict | list | set | BaseModel):
            fields.append((name, _COPY, value))
        else:
            fields.append((name, _FIXED, value))
    t = (fields, cls.model_post_init is not BaseModel.model_post_init)
    _TEMPLATES[cls] = t
    return t


def build(cls: Any, **given: Any) -> Any:
    """A model of `cls` with the fields `given` (trusted: not validated) and the other fields'
    defaults, as ``cls.model_construct(**given)`` makes it."""
    if cls.__private_attributes__:  # pragma: no cover - none of qcdata's models has any
        return cls.model_construct(**given)
    fields, post_init = _TEMPLATES.get(cls) or _template(cls)
    values: dict[str, Any] = {}
    for name, kind, default in fields:
        if name in given:
            values[name] = given[name]
        elif kind == _FIXED:
            values[name] = default
        elif kind == _EMPTY:
            values[name] = default()
        elif kind == _COPY:
            values[name] = copy.deepcopy(default)
    m = cls.__new__(cls)
    object.__setattr__(m, "__dict__", values)
    object.__setattr__(m, "__pydantic_fields_set__", set(given))
    object.__setattr__(m, "__pydantic_extra__", None)
    object.__setattr__(m, "__pydantic_private__", None)
    if post_init:
        m.model_post_init(None)
    return m
