// SPDX-License-Identifier: LGPL-3.0-or-later
//! Writes `python/axtb/_build.py`: whether this build of the extension module assumes an
//! x86-64-v3 CPU (AVX2, FMA, ...; `.cargo/config.toml`), so that `axtb/__init__.py` checks the
//! CPU before it loads the module only when the build needs it (a build with
//! `RUSTFLAGS="-C target-cpu=x86-64"` imports on any x86-64 CPU), and the git commit the build
//! was made from (`provenance.program_commit`; `None` outside a git checkout, e.g. from the
//! source distribution).

use std::path::{Path, PathBuf};
use std::process::Command;

/// The output of `git args` in `dir`, trimmed, if git runs and succeeds.
fn git(dir: &Path, args: &[&str]) -> Option<String> {
    let out = Command::new("git").args(args).current_dir(dir).output().ok()?;
    out.status.success().then(|| String::from_utf8_lossy(&out.stdout).trim().to_string()).filter(|s| !s.is_empty())
}

fn main() {
    let arch = std::env::var("CARGO_CFG_TARGET_ARCH").unwrap_or_default();
    let features = std::env::var("CARGO_CFG_TARGET_FEATURE").unwrap_or_default();
    let v3 = arch == "x86_64" && features.split(',').any(|f| f == "avx2" || f == "fma");
    let dir = PathBuf::from(std::env::var("CARGO_MANIFEST_DIR").expect("cargo sets CARGO_MANIFEST_DIR"));
    let commit = git(&dir, &["rev-parse", "HEAD"]);
    let text = format!(
        "# Written by crates/axtb-py/build.rs (not tracked): whether the extension module was built for\n\
         # x86-64-v3 CPUs (AVX2, FMA, ...), which axtb/__init__.py then checks for before loading it, and the\n\
         # git commit of the build (None outside a git checkout).\n\
         X86_64_V3 = {}\n\
         COMMIT = {}\n",
        if v3 { "True" } else { "False" },
        commit.as_deref().map_or("None".to_string(), |c| format!("\"{c}\""))
    );
    // python/ next to this file in the repository; in maturin's source distribution the python
    // source sits at the top, next to the workspace's Cargo.toml
    let package = ["python/axtb", "../../python/axtb"]
        .iter()
        .map(|p| dir.join(p))
        .find(|p| p.is_dir())
        .expect("python/axtb next to crates/axtb-py/Cargo.toml or two directories up");
    let path = package.join("_build.py");
    // rewritten only when it changes, so that an unchanged build leaves the file alone
    if std::fs::read_to_string(&path).ok().as_deref() != Some(text.as_str()) {
        std::fs::write(&path, text).expect("write python/axtb/_build.py");
    }
    // rerun for a new target CPU (RUSTFLAGS), a deleted file or a new commit, not for every edit
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-changed={}", path.display());
    println!("cargo:rerun-if-env-changed=RUSTFLAGS");
    println!("cargo:rerun-if-env-changed=CARGO_ENCODED_RUSTFLAGS");
    // the HEAD file and the branch it names (when not packed): a new commit reruns this script
    let mut watch = vec![git(&dir, &["rev-parse", "--path-format=absolute", "--git-path", "HEAD"])];
    if let Some(r) = git(&dir, &["symbolic-ref", "-q", "HEAD"]) {
        watch.push(git(&dir, &["rev-parse", "--path-format=absolute", "--git-path", &r]));
    }
    for p in watch.into_iter().flatten().filter(|p| Path::new(p).exists()) {
        println!("cargo:rerun-if-changed={p}");
    }
}
