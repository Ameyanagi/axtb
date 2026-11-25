// SPDX-License-Identifier: LGPL-3.0-or-later

mod geometry;
mod potential;
mod runner;

use std::path::PathBuf;

use anyhow::Result;
use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(about = "Benchmark ReShiki geometry optimization on axtb and COSMolKit")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Generate deterministic starts from a case manifest.
    Prepare {
        #[arg(long)]
        manifest: PathBuf,
        #[arg(long)]
        out: PathBuf,
    },
    /// Optimize selected starts for each configuration.
    Run(runner::RunArgs),
    /// Refine the lowest same-surface candidate at tight thresholds.
    Reference(runner::ReferenceArgs),
}

fn main() -> Result<()> {
    // SAFETY: this is the first operation in main, before starting worker threads
    // or calling any dependency that could read the process environment.
    unsafe {
        for name in [
            "RAYON_NUM_THREADS",
            "OMP_NUM_THREADS",
            "OMP_THREAD_LIMIT",
            "OPENBLAS_NUM_THREADS",
            "MKL_NUM_THREADS",
            "BLIS_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS",
        ] {
            std::env::set_var(name, "1");
        }
    }
    match Cli::parse().command {
        Command::Prepare { manifest, out } => geometry::prepare(&manifest, &out),
        Command::Run(args) => runner::run(args),
        Command::Reference(args) => runner::reference(args),
    }
}
