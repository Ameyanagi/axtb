# SPDX-License-Identifier: LGPL-3.0-or-later
"""The axtb JSON CLI and provenance shared by both validation scripts."""

import hashlib
import json
import os
import platform
import shutil
import subprocess
from pathlib import Path

from cases import ATOMIC_NUMBERS, BOHR_2018, BOHR_2022, ROOT

TIMEOUT = 1800  # Seconds per calculation, as in tests/reference_matrix/runners.py.
METHOD_NAMES = {
    **dict.fromkeys(("gfn2", "gfn2xtb"), "gfn2"),
    "gfnff": "gfnff",
    **dict.fromkeys(("gxtbv11", "gxtb11"), "gxtb-v1.1"),
    **dict.fromkeys(("gxtbv201", "gxtb201", "gxtbv20", "gxtb20"), "gxtb-v2.0.1"),
}
CONVENTION_NAMES = {
    **dict.fromkeys(("axtb", "default"), "axtb"),
    **dict.fromkeys(("xtb", "xtb680"), "xtb"),
    "xtb671": "xtb-6.7.1",
    "tblite": "tblite",
    **dict.fromkeys(("gxtb", "gxtbprogram"), "gxtb"),
}
HESSIAN_NAMES = {
    **dict.fromkeys(("analytic", "analytical"), "analytic"),
    "response": "response",
    **dict.fromkeys(("numerical", "numeric", "fd"), "numerical"),
}


class AxtbError(RuntimeError):
    def __init__(self, kind, message, response=None):
        super().__init__(message)
        self.kind = kind
        self.message = message
        self.response = response


def axtb_env():
    env = os.environ.copy()
    for key in ("AXTB_LAPACK", "AXTB_LAPACK_BLAS"):
        env.pop(key, None)
    return env


def axtb_binary():
    value = os.environ.get("AXTB_BIN", str(ROOT / "target/release/axtb"))
    found = shutil.which(os.path.expanduser(value))
    if not found:
        raise RuntimeError(f"axtb CLI (AXTB_BIN) is missing or not executable: {value}")
    return str(Path(found).resolve())


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def axtb_metadata(binary):
    try:
        proc = subprocess.run(
            [binary, "--version"], capture_output=True, text=True, timeout=30, env=axtb_env(), check=False
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise RuntimeError(f"axtb CLI (AXTB_BIN) --version failed: {exc}") from exc
    version = proc.stdout.strip()
    if proc.returncode or not version.startswith("axtb "):
        raise RuntimeError(f"AXTB_BIN did not identify itself as axtb: {binary}: {proc.stdout}{proc.stderr}")
    return {"version": version, "sha256": sha256(binary), "path": Path(binary).name}


def host_metadata():
    cpu = ""
    if platform.system() == "Linux" and Path("/proc/cpuinfo").is_file():
        for line in Path("/proc/cpuinfo").read_text().splitlines():
            if line.partition(":")[0].strip() == "model name":
                cpu = line.partition(":")[2].strip()
                break
    if platform.system() == "Darwin":
        proc = subprocess.run(
            ["/usr/sbin/sysctl", "-n", "machdep.cpu.brand_string"], capture_output=True, text=True, check=False
        )
        cpu = proc.stdout.strip()
    return {
        "platform": platform.platform(),
        "machine": platform.machine(),
        "cpu": cpu or platform.processor() or "unknown",
    }


def _canonical_name(value, names, field, default):
    # Match crates/driver/src/names.rs: only ASCII letters and digits matter.
    if value is None:
        value = default
    if isinstance(value, str):
        key = "".join(char.lower() for char in value if char.isascii() and char.isalnum())
        if key in names:
            return names[key]
    raise AxtbError("parse", f"Successful CLI response to unrecognized {field}: {value!r}")


def _response_identity(payload, result):
    """Verify the resolved settings before trusting any returned quantities."""
    analysis = payload.get("analyze_hessian") is not None
    method = None if analysis else _canonical_name(payload.get("method"), METHOD_NAMES, "method", "gfn2")
    convention_field = "analysis_conventions" if analysis else "conventions"
    conventions = _canonical_name(payload.get(convention_field), CONVENTION_NAMES, convention_field, "axtb")
    if analysis and conventions == "gxtb":
        conventions = "xtb-6.7.1"
    expected = {"method": method, "conventions": conventions}
    if payload.get("hessian"):
        default = "numerical" if method == "gxtb-v1.1" else "analytic"
        expected["hessian_method"] = _canonical_name(
            payload.get("hessian_method"), HESSIAN_NAMES, "hessian_method", default
        )
    for field, value in expected.items():
        if field not in result or result[field] != value:
            raise AxtbError(
                "parse", f"Expected response {field} {value!r}, got {result.get(field)!r}", {field: result.get(field)}
            )
    return method


def request(binary, payload):
    try:
        proc = subprocess.run(
            [binary, "--json"],
            input=json.dumps(payload, allow_nan=False),
            capture_output=True,
            text=True,
            env=axtb_env(),
            timeout=TIMEOUT,
            check=False,
        )
    except subprocess.TimeoutExpired as exc:
        raise AxtbError("timeout", str(exc)) from exc
    except OSError as exc:
        raise AxtbError("process", str(exc)) from exc
    try:
        result = json.loads(proc.stdout)
    except ValueError as exc:
        raise AxtbError("process" if proc.returncode else "parse", proc.stdout + proc.stderr) from exc
    if not isinstance(result, dict):
        raise AxtbError("parse", proc.stdout + proc.stderr)
    if "error" in result:
        error = result["error"]
        if not isinstance(error, dict) or not {"kind", "message"} <= error.keys():
            raise AxtbError("parse", proc.stdout + proc.stderr)
        raise AxtbError(error["kind"], error["message"], error)
    if proc.returncode:
        raise AxtbError("process", proc.stdout + proc.stderr)
    method = _response_identity(payload, result)
    eigensolver = result.get("eigensolver")
    needs_solver = method is not None and method != "gfnff"
    if eigensolver != "faer" and (needs_solver or eigensolver is not None):
        raise AxtbError("eigensolver", f"Expected faer, got {eigensolver!r}", {"eigensolver": eigensolver})
    return result


def calculation_payload(structure, charge, uhf, method, *, positions_bohr=None, **settings):
    symbols, positions = structure
    if positions_bohr is None:
        bohr = BOHR_2022 if settings.get("conventions", "axtb") == "axtb" else BOHR_2018
        positions_bohr = positions / bohr
    return {
        "numbers": [ATOMIC_NUMBERS[s] for s in symbols],
        "positions": positions_bohr.tolist(),
        "charge": charge,
        "uhf": uhf,
        "method": "gxtb-v2.0.1" if method == "gxtb" else method,
        "threads": 1,
        **settings,
    }


def run_axtb(binary, structure, charge, uhf, method, *, positions_bohr=None, **settings):
    return request(
        binary,
        calculation_payload(structure, charge, uhf, method, positions_bohr=positions_bohr, **settings),
    )
