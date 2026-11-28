#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Run the x86 packing regression with a synthetic 260 MiB L3 cache.

Builds the vendored kernel as a dependency of a temporary test program, so its
BLAS-dependent development tests are not needed. The cache override is confined
to the temporary copy. --kernel-source can select an older copy for comparison.
"""

import argparse
import json
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--kernel-source", type=Path, default=ROOT / "third_party/private-gemm-x86"
    )
    parser.add_argument("--target", help="Optional Cargo target for cross compilation")
    parser.add_argument("--build-only", action="store_true")
    parser.add_argument(
        "--binary", type=Path, help="Copy the built regression executable here"
    )
    options = parser.parse_args()
    if options.build_only and options.binary is None:
        parser.error("--build-only needs --binary")

    with tempfile.TemporaryDirectory(prefix="axtb-packing-budget-") as temporary:
        directory = Path(temporary)
        kernel = directory / "kernel"
        shutil.copytree(
            options.kernel_source,
            kernel,
            ignore=shutil.ignore_patterns("target", "Cargo.lock"),
        )
        cache = kernel / "src/cache.rs"
        text, count = re.subn(
            r"(?ms)^fn cache_info\(\) -> Option<\[CacheInfo; 3\]> \{.*?^\}",
            """fn cache_info() -> Option<[CacheInfo; 3]> {
    Some([
        CacheInfo { associativity: 8, cache_bytes: 32 * 1024, cache_line_bytes: 64 },
        CacheInfo { associativity: 8, cache_bytes: 2 * 1024 * 1024, cache_line_bytes: 64 },
        CacheInfo { associativity: 8, cache_bytes: 260 * 1024 * 1024, cache_line_bytes: 64 },
    ])
}""",
            cache.read_text(),
            count=1,
        )
        if count != 1:
            raise ValueError(
                "Expected exactly one cache_info function in the test copy"
            )
        cache.write_text(text)
        manifest = directory / "Cargo.toml"
        manifest.write_text(
            """[package]
name = "axtb-packing-budget"
version = "0.0.0"
edition = "2021"
[workspace]
exclude = ["kernel"]
[dependencies]
private-gemm-x86 = { path = "kernel", features = ["std", "rayon"] }
rayon = "1"
[[bin]]
name = "packing-budget"
path = "main.rs"
[profile.release]
opt-level = 2
"""
        )
        shutil.copyfile(
            ROOT / "third_party/private-gemm-x86/tests/packing_budget.rs",
            directory / "main.rs",
        )
        target = Path(
            os.environ.get("CARGO_TARGET_DIR", ROOT / "target/packing-budget")
        ).resolve()
        command = ["cargo", "build", "--release", "--manifest-path", str(manifest)]
        if options.target:
            command.extend(["--target", options.target])
        environment = dict(os.environ, CARGO_TARGET_DIR=str(target))
        print(
            "Packing regression: synthetic 260 MiB L3, 256 MiB total heap", flush=True
        )
        subprocess.run(command, env=environment, check=True)
        binary = target
        if options.target:
            binary /= options.target
        binary = (
            binary
            / "release"
            / ("packing-budget.exe" if os.name == "nt" else "packing-budget")
        )
        if options.binary:
            options.binary.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(binary, options.binary)
            options.binary.chmod(0o755)
        if not options.build_only:
            subprocess.run([binary], check=True)
        print(
            json.dumps(
                {"cache_mib": 260, "heap_mib": 256, "built_only": options.build_only}
            )
        )


if __name__ == "__main__":
    main()
