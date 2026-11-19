# SPDX-License-Identifier: LGPL-3.0-or-later
"""axtb's keywords (``axtb._keywords``): the strict checks in their order, the defaults and the
resolved settings (which reproduce the plan when passed back), the driver request built from
the settings given, the 0.1 names' suggestions, the schemas and the generated reference."""

import json
from pathlib import Path

import numpy as np
import pytest

qcdata = pytest.importorskip("qcdata", minversion="0.24.0a5")
from qcdata import ProgramError, ProgramInput, Structure  # noqa: E402

from axtb import _axtb, _keywords as K  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
WATER = Structure(symbols=["O", "H", "H"], geometry=[[0, 0, 0.22], [0, 1.43, -0.89], [0, -1.43, -0.89]])
O2 = Structure(symbols=["O", "O"], geometry=[[0, 0, 0], [0, 0, 2.28]], multiplicity=3)
SCF = ("gfn2", "gxtb-v1.1", "gxtb-v2.0.1")
ALL_CONVENTIONS = [(m, c) for m in K.METHODS for c in K.CONVENTIONS[m]]


def inp(method="gfn2", calctype="energy", keywords=None, structure=WATER, **kw):
    return ProgramInput(program="axtb", calctype=calctype, structure=structure, model={"method": method},
                        keywords=keywords or {}, **kw)


def plan(i, **kw):
    return K.plan_input(i, **kw)


def refused(i, **kw):
    """(kind, code, message) of the refusal of `i`."""
    with pytest.raises(ProgramError) as e:
        K.plan_input(i, **kw)
    return e.value.error.kind.value, e.value.error.code, e.value.error.message


def nested(path, value):
    group, _, leaf = path.rpartition(".")
    return {group: {leaf: value}} if group else {leaf: value}


# ----------------------------------------------------------------------------- names


@pytest.mark.parametrize("name,canonical", [("GFN2-xTB", "gfn2"), ("gfn-ff", "gfnff"), ("g-xTB 2.0", "gxtb-v2.0.1"),
                                            ("gxtb_v1.1", "gxtb-v1.1")])
def test_method_aliases(name, canonical):
    assert plan(inp(name)).method == canonical


@pytest.mark.parametrize("name", ["gxtb", "g-xTB", "gfn3", ""])
def test_unknown_methods(name):
    kind, code, msg = refused(inp(name))
    assert (kind, code) == ("input", "invalid_value")
    if "xtb" in name.lower():
        assert "gxtb-v2.0.1" in msg and "gxtb-v1.1" in msg


@pytest.mark.parametrize("method,given,canonical", [
    ("gfn2", "XTB", "xtb-6.8.0"), ("gfn2", "xtb 6.7.1", "xtb-6.7.1"), ("gfn2", "tblite", "tblite-0.8.0"),
    ("gxtb-v1.1", "gxtb", "gxtb-1.1.0"), ("gxtb-v2.0.1", "gxtb", "gxtb-2.0.1"), ("gxtb-v2.0.1", "G-XTB 2.0.1",
                                                                                   "gxtb-2.0.1"),
    ("gfnff", "xtb", "xtb-6.8.0"), ("gfnff", "AXTB", "axtb"),
])
def test_conventions_names(method, given, canonical):
    p = plan(inp(method, keywords={"conventions": given}))
    assert p.settings["conventions"] == canonical
    assert p.request["conventions"] == K.CONVENTIONS[method][canonical]


@pytest.mark.parametrize("method,given", [("gxtb-v1.1", "gxtb-2.0.1"), ("gfnff", "tblite"), ("gfn2", "default")])
def test_conventions_of_another_method(method, given):
    assert refused(inp(method, keywords={"conventions": given}))[:2] == ("input", "invalid_value")


def test_basis_is_refused():
    i = ProgramInput(program="axtb", calctype="energy", structure=WATER, model={"method": "gfn2", "basis": "sto-3g"})
    kind, code, msg = refused(i)
    assert (kind, code) == ("input", "invalid_value") and msg.startswith("model.basis")


# ----------------------------------------------------------------------------- strictness


@pytest.mark.parametrize("old", sorted(K.RENAMED))
def test_every_old_name_suggests_the_new_one(old):
    kind, code, msg = refused(inp("gfn2", "hessian", {old: 1}))
    if old in K.GROUPS:  # a group given a value
        assert (kind, code) == ("input", "invalid_value") and "0.1 key" in msg
        return
    assert (kind, code) == ("input", "unknown_key")
    assert "did you mean" in msg
    target = K.RENAMED[old].split(" ")[0].rstrip(",")
    assert target.split("(")[0] in msg


def test_unknown_keys():
    assert refused(inp(keywords={"zz": 1}))[:2] == ("input", "unknown_key")
    kind, code, msg = refused(inp(keywords={"scf": {"maxiter": 3}}))
    assert (kind, code) == ("input", "unknown_key") and "scf.max_iter" in msg
    assert "scf.etemp" in refused(inp(keywords={"scf": {"etmp": 300}}))[2]
    # a group given a scalar
    assert refused(inp(keywords={"scf": 3}))[:2] == ("input", "invalid_value")


def test_order_unknown_before_invalid_before_no_effect():
    # an unknown key wins over an invalid value
    assert refused(inp(keywords={"scf": {"max_iter": 0}, "zz": 1}))[1] == "unknown_key"
    # an invalid value wins over no_effect
    assert refused(inp("gfnff", keywords={"scf": {"max_iter": 0}}))[1] == "invalid_value"
    # no_effect wins over not_applicable
    assert refused(inp("gfnff", "energy", {"scf": {"max_iter": 10}, "hessian": {"method": "analytic"}}))[1] == "no_effect"
    # not_applicable wins over missing_setting
    assert refused(inp("gfn2", "energy", {"solvent": {"model": "gbsa"}, "freq": {"pressure": 1e5}}))[1] == "not_applicable"
    # missing_setting wins over the structure
    bad = Structure(symbols=["O", "H"], geometry=[[0, 0, 0], [0, 0, 1.8]])
    assert refused(inp("gfn2", "energy", {"solvent": {"model": "gbsa"}}, structure=bad))[1] == "missing_setting"
    # the structure wins over unsupported settings
    assert refused(inp("gfn2", "energy", {"scf": {"e_tol": 1e-6}}, structure=bad))[1] == "invalid_structure"


@pytest.mark.parametrize("keywords,code", [
    ({"threads": 0}, "invalid_value"),
    ({"threads": True}, "invalid_value"),
    ({"threads": 1.5}, "invalid_value"),
    ({"accuracy_scale": -1}, "invalid_value"),
    ({"accuracy": 0.1}, "invalid_value"),
    ({"accuracy": {"gradient": -1}}, "invalid_value"),
    ({"accuracy": {"grad": 1e-5}}, "unknown_key"),
    ({"accuracy": {}}, "no_effect"),
    ({"scf": {"max_iter": 0}}, "invalid_value"),
    ({"scf": {"etemp": -1}}, "invalid_value"),
    ({"scf": {"etemp": float("nan")}}, "invalid_value"),
    ({"scf": {"mixer": "anderson"}}, "invalid_value"),
    ({"scf": {"guess": "eeq-bc"}}, "invalid_value"),
    ({"scf": {"restart": "yes"}}, "invalid_value"),
    ({"solvent": {"name": "watre"}}, "invalid_value"),
    ({"solvent": {"name": "water", "model": "pcm"}}, "invalid_value"),
    ({"solvent": {"name": "water", "state": "1M"}}, "invalid_value"),
    ({"output": {"orbitals": 1}}, "invalid_value"),
])
def test_invalid_values(keywords, code):
    assert refused(inp("gfn2", "energy", keywords))[:2] == ("input", code)


def test_solvent_did_you_mean():
    msg = refused(inp(keywords={"solvent": {"name": "watre"}}))[2]
    assert "water" in msg


@pytest.mark.parametrize("method,keywords,path", [
    ("gfnff", {"accuracy_scale": 0.1}, "accuracy_scale"),
    ("gfnff", {"scf": {"etemp": 300}}, "scf.etemp"),
    ("gfnff", {"output": {"orbitals": False}}, "output.orbitals"),
    ("gfn2", {"topology": {"fragment_charges": [0]}}, "topology.fragment_charges"),
    ("gxtb-v2.0.1", {"scf": {"mixer": "diis"}}, "scf.mixer"),
    ("gxtb-v2.0.1", {"solvent": {"name": "water", "correction": "gfn2"}, "scf": {"guess": "eeq-bc"}}, "scf.guess"),
    ("gfn2", {"solvent": {"name": "water", "correction": "gfn2"}}, "solvent.correction"),
    ("gxtb-v2.0.1", {"conventions": "gxtb", "output": {"hamiltonian": False}}, "output.hamiltonian"),
    ("gfn2", {"conventions": "tblite", "scf": {"guess": "eeq"}}, "scf.guess"),
    ("gfn2", {"scf": {"restart": "off"}}, "scf.restart"),
])
def test_no_effect(method, keywords, path):
    kind, code, msg = refused(inp(method, "energy", keywords))
    assert (kind, code) == ("input", "no_effect") and msg.startswith(path)


@pytest.mark.parametrize("keywords,path", [
    ({"hessian": {"method": "analytic", "step": 0.01}}, "hessian.step"),
    ({"hessian": {"step": 0.01}}, "hessian.step"),
    ({"freq": {"analysis": "off", "temperature": 300}}, "freq.temperature"),
    ({"freq": {"thermo": "off", "pressure": 1e5}}, "freq.pressure"),
])
def test_no_effect_of_other_settings(keywords, path):
    kind, code, msg = refused(inp("gfn2", "hessian", keywords))
    assert (kind, code) == ("input", "no_effect") and msg.startswith(path)


@pytest.mark.parametrize("calctype", ["energy", "gradient"])
@pytest.mark.parametrize("path,value", [("hessian.method", "numerical"), ("freq.temperature", 300.0),
                                        ("freq.analysis", "off")])
def test_hessian_settings_need_a_hessian(calctype, path, value):
    kind, code, msg = refused(inp("gfn2", calctype, nested(path, value)))
    assert (kind, code) == ("input", "not_applicable") and msg.startswith(path)


def test_session_calctypes_decide_applicability():
    kw = {"hessian": {"method": "numerical"}}
    assert refused(inp("gfn2", "gradient", kw), calctypes=frozenset({"energy", "gradient"}), session=True)[1] == \
        "not_applicable"
    p = plan(inp("gfn2", "gradient", kw), calctypes=frozenset({"energy", "gradient", "hessian"}), session=True)
    assert "hessian.method" in p.contingent and p.task["hessian_method"] == "numerical"


@pytest.mark.parametrize("field,value", [("cmdline_args", ["--x"]), ("files", {"a": "b"})])
def test_input_fields_axtb_takes_none_of(field, value):
    assert refused(inp(**{field: value}))[:2] == ("input", "not_applicable")


@pytest.mark.parametrize("method,keywords,path", [
    ("gfn2", {"solvent": {"model": "gbsa"}}, "solvent.model"),
    ("gfn2", {"solvent": {"state": "bar1M"}}, "solvent.state"),
    ("gxtb-v2.0.1", {"solvent": {"correction": "gfn2"}}, "solvent.correction"),
    ("gxtb-v2.0.1", {"solvent": {"name": "water"}}, "solvent.name"),
    ("gxtb-v1.1", {"solvent": {"correction_conventions": "axtb"}}, "solvent.correction_conventions"),
    ("gfn2", {"scf": {"guess": "given"}}, "scf.guess"),
])
def test_missing_settings(method, keywords, path):
    kind, code, msg = refused(inp(method, "energy", keywords))
    assert (kind, code) == ("input", "missing_setting") and msg.startswith(path)


def test_gxtb_solvent_hint():
    msg = refused(inp("gxtb-v2.0.1", keywords={"solvent": {"name": "water"}}))[2]
    assert 'solvent.correction "gfn2"' in msg and "composite" in msg


def test_structures():
    radical = Structure(symbols=["O", "H"], geometry=[[0, 0, 0], [0, 0, 1.8]])
    kind, code, msg = refused(inp("gfn2", structure=radical))
    assert (kind, code) == ("input", "invalid_structure") and "2" in msg
    plan(inp("gfn2", structure=radical.model_copy(update={"multiplicity": 2})))
    assert refused(inp("gfnff", structure=O2))[:2] == ("input", "not_modeled")
    plan(inp("gfnff", structure=O2.model_copy(update={"multiplicity": 1})))
    nd = Structure(symbols=["Nd", "O"], geometry=[[0, 0, 0], [0, 0, 3.5]], multiplicity=1)
    with pytest.raises(ProgramError) as e:
        plan(inp("gxtb-v1.1", structure=nd))
    assert (e.value.error.kind.value, e.value.error.code) == ("input", "unsupported_element")
    # the composite needs GFN2-xTB's elements too (H-Rn), g-xTB 2.0.1 alone has H-Lr
    fr = Structure(symbols=["Fr", "H"], geometry=[[0, 0, 0], [0, 0, 4.0]])
    plan(inp("gxtb-v2.0.1", structure=fr))
    assert refused(inp("gxtb-v2.0.1", keywords={"solvent": {"name": "water", "correction": "gfn2"}},
                       structure=fr))[1] == "unsupported_element"


@pytest.mark.parametrize("keywords,code", [
    ({"scf": {"e_tol": 1e-8}}, "unsupported_setting"),
    ({"scf": {"g_tol": 1e-8}}, "unsupported_setting"),
    ({"solvent": {"name": "water", "eps": 80}}, "unsupported_setting"),
    ({"output": {"molden": "spherical"}}, "unsupported_setting"),
    ({"output": {"density": True}}, "unsupported_setting"),
    ({"output": {"core_hamiltonian": True}}, "not_applicable"),
    ({"memory_gib": 4}, "unsupported_setting"),
    ({"backend": "gpu"}, "unsupported_setting"),
    ({"device_memory_gib": 4}, "unsupported_setting"),
])
def test_unsupported(keywords, code):
    assert refused(inp("gfn2", "energy", keywords))[:2] == ("unsupported", code)


@pytest.mark.parametrize("keywords", [{"output": {"orbitals": True}}, {"output": {"bond_orders": True}},
                                      {"output": {"hamiltonian": True}}, {"solvent": {"name": "water"}}])
def test_gfnff_refuses_what_it_cannot_give(keywords):
    assert refused(inp("gfnff", "energy", keywords))[0] == "unsupported"


def test_hamiltonian_in_program_conventions():
    kind, code, _ = refused(inp("gxtb-v2.0.1", keywords={"conventions": "gxtb", "output": {"hamiltonian": True}}))
    assert (kind, code) == ("unsupported", "not_applicable")
    plan(inp("gxtb-v1.1", keywords={"conventions": "gxtb", "output": {"hamiltonian": True}}))


@pytest.mark.parametrize("calctype", ["optimization", "transition_state", "irc", "conformer_search"])
def test_unsupported_calctypes(calctype):
    kind, code, msg = refused(inp(calctype=calctype))
    assert (kind, code) == ("unsupported", "unsupported_calctype") and "qcopt" in msg


@pytest.mark.parametrize("method,conventions,hm,ok", [
    ("gxtb-v1.1", "axtb", "analytic", False), ("gxtb-v1.1", "axtb", "response", True),
    ("gxtb-v1.1", "gxtb", "response", False), ("gfnff", "axtb", "response", False),
    ("gfnff", "axtb", "numerical", True), ("gxtb-v2.0.1", "gxtb", "response", True),
])
def test_hessian_methods_offered(method, conventions, hm, ok):
    i = inp(method, "hessian", {"conventions": conventions, "hessian": {"method": hm}})
    if ok:
        assert plan(i).task["hessian_method"] == hm
    else:
        assert refused(i)[:2] == ("input", "invalid_value")


# ----------------------------------------------------------------------------- defaults and resolved


@pytest.mark.parametrize("method,conventions", ALL_CONVENTIONS)
@pytest.mark.parametrize("calctype", ["energy", "hessian"])
def test_resolved_settings_reproduce_the_plan(method, conventions, calctype):
    p = plan(inp(method, calctype, {"conventions": conventions}))
    again = plan(inp(method, calctype, p.settings))
    assert again.settings == p.settings


@pytest.mark.parametrize("session", [False, True])
@pytest.mark.parametrize("keywords", [
    {"solvent": {"name": "water", "correction": "gfn2"}},
    {"conventions": "gxtb", "solvent": {"name": "toluene", "model": "gbsa", "correction": "gfn2"}},
    {"conventions": "gxtb", "solvent": {"name": "water", "correction": "gfn2", "correction_conventions": "axtb"}},
    {"scf": {"etemp": 1000}, "solvent": {"name": "water", "correction": "gfn2"}},
    {"scf": {"etemp": 1000, "retries": [{"etemp": 3000}]}, "solvent": {"name": "water", "correction": "gfn2"}},
])
@pytest.mark.parametrize("method", ["gxtb-v1.1", "gxtb-v2.0.1"])
def test_composite_resolved(method, keywords, session):
    p = plan(inp(method, "hessian", keywords), session=session)
    assert p.composite
    assert plan(inp(method, "hessian", p.settings), session=session).settings == p.settings
    if "etemp" not in json.dumps(keywords):
        assert "etemp" not in p.settings["scf"] and "scf.etemp" in p.reasons
    else:
        assert p.request["etemp"] == 1000.0
    if method == "gxtb-v1.1":
        assert "max_iter" not in p.settings["scf"]
    assert p.request["solvation_correction"] == "gfn2"


@pytest.mark.parametrize("method,etemp,max_iter,mixer,guess,hm,step,hacc", [
    ("gfn2", 300.0, 250, "broyden-then-diis", "eeq", "analytic", None, 1e-3),
    ("gxtb-v2.0.1", 0.0, 250, None, "eeq-bc", "analytic", None, 1e-5),
    ("gxtb-v1.1", 0.0, 1000, None, "eeq-bc", "numerical", 0.005, 0.01),
])
def test_defaults(method, etemp, max_iter, mixer, guess, hm, step, hacc):
    s = plan(inp(method, "hessian")).settings
    assert s["scf"]["etemp"] == etemp and s["scf"]["max_iter"] == max_iter and s["scf"]["guess"] == guess
    assert s["scf"].get("mixer") == mixer
    assert s["hessian"]["method"] == hm and s["hessian"].get("step") == step
    assert s["hessian"]["accuracy_scale"] == hacc
    assert s["freq"] == {"analysis": "harmonic", "conventions": "axtb", "thermo": "qrrho", "temperature": 298.15,
                         "pressure": 101325.0, "scale": 1.0, "imaginary_cutoff": 20.0, "rotor_cutoff": 50.0,
                         "symmetry_tol": 0.1}


def test_program_conventions_defaults():
    s = plan(inp("gxtb-v2.0.1", "hessian", {"conventions": "gxtb"}), session=True).settings
    assert s["freq"]["conventions"] == "xtb-6.7.1" and s["scf"]["restart"] == "off"
    assert "hamiltonian" not in s["output"]
    s = plan(inp("gfnff", "hessian", {"conventions": "xtb"})).settings
    assert s["freq"]["scale"] == 1.03 and s["freq"]["conventions"] == "xtb-6.8.0"
    s = plan(inp("gfn2", "energy", {"conventions": "tblite"})).settings
    assert s["scf"]["guess"] == "sad" and s["scf"]["mixer"] == "broyden"


def test_driver_request_holds_only_what_was_given():
    p = plan(inp("gfn2", "hessian"))
    assert p.request == {"numbers": [], "positions": [], "method": "gfn2", "charge": 0.0, "uhf": 0}
    assert p.task == {}
    kw = {"scf": {"etemp": 1000, "max_iter": 77, "mixer": "broyden-then-diis", "guess": "sad"},
          "solvent": {"name": "h2o", "model": "gbsa", "state": "bar1m"}, "accuracy_scale": 0.1,
          "hessian": {"method": "numerical", "step": 0.01, "accuracy_scale": 0.001},
          "freq": {"conventions": "xtb", "temperature": [200, 300], "symmetry_number": 2}}
    p = plan(inp("gfn2", "hessian", kw, structure=O2))
    assert p.request == {"numbers": [], "positions": [], "method": "gfn2", "charge": 0.0, "uhf": 2,
                         "etemp": 1000.0, "max_iter": 77, "mixer": "auto", "guess": "sad", "gbsa": "water",
                         "solvation_state": "bar1M", "accuracy": 0.1}
    assert p.task == {"hessian_method": "numerical", "hessian_step": 0.01, "hessian_accuracy": 0.001,
                      "analysis_conventions": "xtb", "temperature": [200.0, 300.0], "symmetry_number": 2}
    assert p.settings["solvent"] == {"name": "water", "model": "gbsa", "state": "bar1M"}
    # a session keeps the restart policy beside a named start
    p = plan(inp("gfn2", "gradient", {"scf": {"guess": "sad"}}), session=True)
    assert (p.request["guess"], p.request["restart"]) == ("sad", "previous")
    p = plan(inp("gxtb-v2.0.1", "gradient", {"scf": {"restart": "previous"}}), session=True)
    assert p.request["restart"] == "previous"


@pytest.mark.parametrize("method", SCF)
def test_retries_resolved_and_requested(method):
    kw = {"scf": {"retries": [{"etemp": 1000, "max_iter": 1000}, {"max_iter": 500}]}}
    p = plan(inp(method, "gradient", kw))
    assert p.settings["scf"]["retries"] == [{"etemp": 1000.0, "max_iter": 1000}, {"max_iter": 500}]
    assert p.request["scf_retries"] == p.settings["scf"]["retries"]
    assert "scf.retries" in p.record["contingent"]
    assert plan(inp(method, "gradient", p.settings)).settings == p.settings
    # a reused plan's rungs are copies
    p.settings["scf"]["retries"][0]["etemp"] = 5.0
    assert plan(inp(method, "gradient", kw)).settings["scf"]["retries"][0]["etemp"] == 1000.0
    assert "retries" not in plan(inp(method, "gradient")).settings["scf"]


def test_gfnff_topology_and_accuracy():
    geo = (np.asarray(WATER.geometry) + 0.01).tolist()
    p = plan(inp("gfnff", "gradient", {"topology": {"geometry": geo, "fragment_charges": [0]},
                                        "accuracy": {"gradient": 1e-6}}))
    assert p.request["topology_positions"] == geo and p.request["fragment_charges"] == [0.0]
    assert p.accuracy.gradient == 1e-6 and "uhf" not in p.request
    assert refused(inp("gfnff", keywords={"topology": {"geometry": [[0, 0, 0]]}}))[1] == "invalid_value"
    assert refused(inp("gfnff", keywords={"topology": {"fragment_charges": [1]}}))[1] == "invalid_value"


def test_thread_cap():
    """axtb's internal cap (axtb_driver::THREAD_CAPS): energies and gradients by size, at most
    8; Hessians uncapped above the one-thread size; the plan has no resources.threads and records
    the cap."""
    big = Structure(symbols=["H", "H"] * 160, geometry=[[0, 0, 3.0 * i + 1.4 * (j % 2)] for i in range(160)
                                                        for j in range(2)])
    assert K.thread_cap(plan(inp("gfn2", "gradient")), "gradient") == 1  # water
    assert K.thread_cap(plan(inp("gfn2", "gradient", structure=big)), "gradient") == 8  # 320 atoms
    assert K.thread_cap(plan(inp("gfn2", "hessian", structure=big)), "hessian") is None
    assert K.thread_cap(plan(inp("gfnff", "gradient", structure=big)), "gradient") == 1
    assert K.thread_cap(plan(inp("gfn2", "hessian")), "hessian") == 1  # water: below the one-thread size
    assert K.thread_cap(plan(inp("gxtb-v2.0.1", "hessian")), "hessian") is None
    assert "| g-xTB 2.0.1 | 1 (below 60), 4 (60-299), 8 (from 300) |" in K.thread_caps()
    from axtb import program

    out = program.plan(inp("gfn2", "gradient", {"threads": 6}))
    assert "resources" not in out and out["record"]["thread_cap"] == 1
    assert out["resolved"]["settings"]["threads"] == 6


# ----------------------------------------------------------------------------- schema and docs


@pytest.mark.parametrize("calctype", K.CALCTYPES)
def test_schema_validates_the_resolved_settings(calctype):
    jsonschema = pytest.importorskip("jsonschema")
    schema = K.settings_schema(calctype)
    jsonschema.validators.validator_for(schema).check_schema(schema)
    validator = jsonschema.Draft202012Validator(schema)
    for method, conventions in ALL_CONVENTIONS:
        s = plan(inp(method, calctype, {"conventions": conventions})).settings
        errors = [e.message for e in validator.iter_errors(s)]
        assert not errors, (method, conventions, errors)
    for key in K.KEYS:
        group, _, leaf = key.path.rpartition(".")
        props = schema["properties"].get(group, {}).get("properties", {}) if group else schema["properties"]
        if key.unsupported is not None or calctype not in key.calctypes:
            assert leaf not in props, key.path
        else:
            assert {"x-tier", "x-concern", "x-applies-to"} <= set(props[leaf]), key.path


def test_generated_reference_is_current():
    text = (ROOT / "docs/python.md").read_text()
    begin, end = "<!-- keywords: generated by scripts/gen_python_docs.py -->", "<!-- end keywords -->"
    assert begin in text and end in text
    assert text.split(begin)[1].split(end)[0].strip() == K.markdown_table().strip(), \
        "docs/python.md is stale: run scripts/gen_python_docs.py"


# the driver's labels of conventions
LABELS = {"axtb": "axtb", "xtb": "xtb-6.8.0", "xtb-6.7.1": "xtb-6.7.1", "tblite": "tblite-0.8.0"}


@pytest.mark.parametrize("method,conventions", ALL_CONVENTIONS)
def test_defaults_are_the_drivers(method, conventions):
    """The table's defaults against what the driver reports for a request without them."""
    p = plan(inp(method, "hessian", {"conventions": conventions}))
    req = json.dumps(p.request)
    task = json.dumps({"numbers": [], "positions": [], "hessian": True, **p.task})
    r = _axtb.compute_once(req, p.numbers, p.positions, task)
    s = p.settings
    assert r["hessian"]["hessian_method"] == s["hessian"]["method"]
    v = r["vibrations"]
    assert LABELS[v["vibrational_conventions"]] == s["freq"]["conventions"]
    t = v["thermo"]
    for key, field in (("pressure", "pressure"), ("rotor_cutoff", "rotor_cutoff"),
                       ("imaginary_cutoff", "imaginary_cutoff"), ("scale", "frequency_scale"),
                       ("temperature", "temperature")):
        assert t[field] == s["freq"][key], key
    if method != "gfnff":
        assert r["scf"]["electronic_temperature"] == s["scf"]["etemp"]
    if "mixer" in s.get("scf", {}):
        assert r["scf"]["mixer_used"] in {"broyden-then-diis": ("broyden", "diis-fallback"),
                                          "broyden": ("broyden",)}[s["scf"]["mixer"]]


@pytest.mark.parametrize("method,conventions", [("gxtb-v2.0.1", "axtb"), ("gxtb-v2.0.1", "gxtb"),
                                                ("gxtb-v1.1", "axtb")])
def test_composite_defaults_are_the_drivers(method, conventions):
    p = plan(inp(method, "gradient", {"conventions": conventions, "solvent": {"name": "water", "correction": "gfn2"}}))
    task = json.dumps({"numbers": [], "positions": [], "gradient": True})
    parts = _axtb.compute_once(json.dumps(p.request), p.numbers, p.positions, task)["composite"]["parts"]
    assert parts["gxtb"]["etemp"] == 0.0 and parts["gfn2_gas"]["etemp"] == 300.0
    if "max_iter" in p.settings["scf"]:
        assert {part["max_iter"] for part in parts.values()} == {p.settings["scf"]["max_iter"]}
    else:
        assert parts["gxtb"]["max_iter"] == 1000 and parts["gfn2_gas"]["max_iter"] == 250
    assert LABELS[parts["gfn2_gas"]["conventions"]] == p.settings["solvent"]["correction_conventions"]


def test_solvent_tables_come_from_the_driver():
    assert "water" in _axtb.solvents("alpb") and "water" in _axtb.solvents("gbsa")
    for name in ("H2O", "Water", "dmso", "acetonitrile"):
        assert plan(inp(keywords={"solvent": {"name": name}})).settings["solvent"]["name"] in _axtb.solvents("alpb")
