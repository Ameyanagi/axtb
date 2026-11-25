// SPDX-License-Identifier: LGPL-3.0-or-later

use std::cell::Cell;
use std::collections::{BTreeMap, HashSet};
use std::fs::{self, File, OpenOptions};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::ops::ControlFlow;
use std::path::{Path, PathBuf};
use std::sync::Mutex;
use std::time::Instant;

use anyhow::{Context, Result, bail, ensure};
use rayon::prelude::*;
#[cfg(feature = "knobs")]
use reshiki_optimize::{Acceleration, InitialHessian, Rebuild, Trust, Update};
use reshiki_optimize::{Accuracy, Level, Method, Options, Potential, Status, optimize};
use serde::{Deserialize, Serialize};

use crate::geometry::{Geometry, Identity, Selection, Start, parent, starts};
use crate::potential::{BenchPotential, Hessian, ScfOptions, Surface, accuracy_name};

struct Config {
    label: String,
    options: Options,
    hessian: Hessian,
    preopt: bool,
}

fn parse_config(label: &str) -> Result<Config> {
    let mut parts = label.split('+');
    let method = match parts.next().unwrap_or_default() {
        "anc" => Method::AncRfo,
        "lbfgs" => Method::Lbfgs,
        #[cfg(feature = "knobs")]
        "internal" => Method::Internal,
        #[cfg(not(feature = "knobs"))]
        "internal" => bail!("internal requires an optimizer with Method::Internal and --features knobs"),
        method => bail!("unknown optimizer method {method:?}"),
    };
    let mut options = Options { method, ..Options::default() };
    let mut hessian = Hessian::Model;
    let mut preopt = false;
    let mut seen = HashSet::new();
    for part in parts {
        let key = part.split_once('=').map_or(part, |(key, _)| key);
        ensure!(seen.insert(key), "repeated config modifier {key:?} in {label:?}");
        match part {
            "tight" => options.level = Level::Tight,
            "notightfinal" => options.tight_final = false,
            "preopt=gfnff" => preopt = true,
            part if part.starts_with("maxsteps=") => {
                options.max_steps = part["maxsteps=".len()..].parse().context("expected a nonnegative maxsteps=<N>")?;
            }
            "hess=model" => {
                hessian = Hessian::Model;
                #[cfg(feature = "knobs")]
                {
                    options.initial_hessian = InitialHessian::Model;
                }
            }
            #[cfg(feature = "knobs")]
            "hess=self" | "hess=gfnff" => {
                options.initial_hessian = InitialHessian::Potential;
                hessian = if part == "hess=self" { Hessian::Same } else { Hessian::Gfnff };
            }
            #[cfg(feature = "knobs")]
            part if part.starts_with("hess=auto") => {
                let max_atoms: usize =
                    part["hess=auto".len()..].parse().context("expected hess=auto<N> with a positive atom limit")?;
                ensure!(max_atoms > 0, "hess=auto<N> needs a positive atom limit");
                options.initial_hessian = InitialHessian::Auto { max_atoms };
                hessian = Hessian::Auto;
            }
            #[cfg(feature = "knobs")]
            "update=bfgs" => options.update = Update::Bfgs,
            #[cfg(feature = "knobs")]
            "update=damped" => options.update = Update::DampedBfgs,
            #[cfg(feature = "knobs")]
            "update=bofill" => options.update = Update::Bofill,
            #[cfg(feature = "knobs")]
            "accel=none" => options.acceleration = Acceleration::None,
            #[cfg(feature = "knobs")]
            "accel=gdiis" => options.acceleration = Acceleration::Gdiis,
            #[cfg(feature = "knobs")]
            "accel=gediis" => options.acceleration = Acceleration::Gediis,
            #[cfg(feature = "knobs")]
            "adaptive" => options.adaptive_accuracy = true,
            #[cfg(feature = "knobs")]
            "trust=global" => options.trust = Trust::Global,
            #[cfg(feature = "knobs")]
            "trust=peratom" => options.trust = Trust::per_atom(),
            #[cfg(feature = "knobs")]
            part if part.starts_with("trust=peratom:") => {
                let (start, max) =
                    part["trust=peratom:".len()..].split_once(':').context("expected trust=peratom:<start>:<max>")?;
                let start: f64 = start.parse().context("invalid start radius in trust=peratom:<start>:<max>")?;
                let max: f64 = max.parse().context("invalid maximum radius in trust=peratom:<start>:<max>")?;
                let Trust::PerAtom { min, atom_cap, .. } = Trust::per_atom() else {
                    unreachable!("Trust::per_atom() must return Trust::PerAtom");
                };
                ensure!(
                    min <= start && start <= max && max.is_finite(),
                    "trust=peratom:<start>:<max> requires {min} <= start <= max with finite radii"
                );
                options.trust = Trust::PerAtom { start, min, max, atom_cap };
            }
            #[cfg(feature = "knobs")]
            part if part.starts_with("reject=") => {
                let ratio: f64 = part["reject=".len()..].parse().context("expected reject=<ratio>")?;
                ensure!(ratio.is_finite() && ratio >= 0., "reject=<ratio> requires a finite nonnegative ratio");
                options.reject_ratio = ratio;
            }
            #[cfg(feature = "knobs")]
            part if part.starts_with("floor=") => {
                let floor: f64 = part["floor=".len()..].parse().context("expected floor=<curvature>")?;
                ensure!(floor.is_finite() && floor > 0., "floor=<curvature> requires finite positive curvature");
                options.curvature_floor = floor;
            }
            #[cfg(feature = "knobs")]
            "rebuild=poor" => options.rebuild = Rebuild::OnPoorPrediction,
            #[cfg(feature = "knobs")]
            part if part.starts_with("rebuild=every") => {
                let interval = part["rebuild=every".len()..].parse().context("expected rebuild=every<N>")?;
                ensure!(interval > 0, "rebuild=every<N> requires a positive interval");
                options.rebuild = Rebuild::Every(interval);
            }
            #[cfg(feature = "knobs")]
            "contacts" => options.contacts = true,
            #[cfg(not(feature = "knobs"))]
            part if matches!(
                key,
                "hess" | "update" | "accel" | "adaptive" | "trust" | "reject" | "floor" | "rebuild" | "contacts"
            ) =>
            {
                bail!("{part:?} requires the extended optimizer API and --features knobs");
            }
            _ => bail!("unknown config modifier {part:?} in {label:?}"),
        }
    }
    Ok(Config { label: label.into(), options, hessian, preopt })
}

#[derive(clap::Args)]
pub struct RunArgs {
    #[arg(long)]
    starts: PathBuf,
    #[arg(long)]
    potential: Surface,
    #[arg(long, value_delimiter = ',', required = true)]
    configs: Vec<String>,
    #[arg(long, default_value = "1", value_parser = clap::value_parser!(u32).range(1..))]
    jobs: u32,
    #[arg(long)]
    out: PathBuf,
    #[command(flatten)]
    selection: Selection,
    #[command(flatten)]
    scf: ScfOptions,
}

#[derive(Clone, Copy, Serialize)]
/// [evaluations, seconds, energy, max gradient, accuracy, gradient norm over free components]
struct Sample(usize, f64, f64, f64, &'static str, f64);

#[derive(Serialize)]
struct Stage {
    potential: Surface,
    status: &'static str,
    error: Option<String>,
    method: Option<String>,
    steps: usize,
    evaluations: usize,
    hessian_evaluations: usize,
    potential_seconds: f64,
    setup_seconds: f64,
    gradient_seconds: f64,
    hessian_seconds: f64,
    optimizer_seconds: f64,
    wall_seconds: f64,
    scf_iterations_total: Option<usize>,
    energy: Option<f64>,
    max_gradient: Option<f64>,
    median_gradient_seconds: Option<f64>,
    #[serde(skip)]
    gradient_times: Vec<f64>,
    #[serde(skip)]
    positions: Option<Vec<f64>>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    trace: Vec<Sample>,
}

impl Stage {
    fn new(potential: Surface) -> Self {
        Self {
            potential,
            status: "error",
            error: None,
            method: None,
            steps: 0,
            evaluations: 0,
            hessian_evaluations: 0,
            potential_seconds: 0.0,
            setup_seconds: 0.0,
            gradient_seconds: 0.0,
            hessian_seconds: 0.0,
            optimizer_seconds: 0.0,
            wall_seconds: 0.0,
            scf_iterations_total: (potential == Surface::Gfn2).then_some(0),
            energy: None,
            max_gradient: None,
            median_gradient_seconds: None,
            gradient_times: Vec::new(),
            positions: None,
            trace: Vec::new(),
        }
    }
}

#[derive(Serialize)]
struct Row<'a> {
    #[serde(flatten)]
    id: &'a Identity,
    potential: Surface,
    config: &'a str,
    atoms: usize,
    status: &'static str,
    error: Option<String>,
    steps: usize,
    evaluations: usize,
    hessian_evaluations: usize,
    potential_seconds: f64,
    setup_seconds: f64,
    gradient_seconds: f64,
    gradient_evaluations: usize,
    hessian_seconds: f64,
    hessian_gradient_equivalents: Option<f64>,
    optimizer_seconds: f64,
    wall_seconds: f64,
    scf_iterations_total: Option<usize>,
    energy: Option<f64>,
    max_gradient: Option<f64>,
    final_xyz: Option<PathBuf>,
    trace: Vec<Sample>,
    preopt: Option<Stage>,
    target: Option<Stage>,
    preopt_gfn2_gradient_equivalents: Option<f64>,
    total_gfn2_gradient_equivalents: Option<f64>,
}

pub fn run(args: RunArgs) -> Result<()> {
    let configs = args.configs.iter().map(|label| parse_config(label)).collect::<Result<Vec<_>>>()?;
    ensure!(args.configs.iter().collect::<HashSet<_>>().len() == configs.len(), "duplicate config labels");
    let starts = starts(&args.starts, &args.selection)?;
    let pool = rayon::ThreadPoolBuilder::new().num_threads(args.jobs as usize).build()?;
    let root = parent(&args.out);
    fs::create_dir_all(root)?;
    let file = OpenOptions::new().create(true).append(true).open(&args.out)?;
    let writer = Mutex::new(BufWriter::new(file));
    pool.install(|| {
        starts.par_iter().flat_map(|start| configs.par_iter().map(move |config| (start, config))).try_for_each(
            |(start, config)| -> Result<()> {
                let row = run_one(start, args.potential, config, args.scf, root)?;
                let line = serde_json::to_string(&row)?;
                let mut writer = writer.lock().map_err(|_| anyhow::anyhow!("result writer was poisoned"))?;
                writeln!(writer, "{line}")?;
                writer.flush()?;
                Ok(())
            },
        )
    })?;
    println!("Wrote {} runs to {}", starts.len() * configs.len(), args.out.display());
    Ok(())
}

fn run_one<'a>(
    start: &'a Start,
    potential: Surface,
    config: &'a Config,
    scf: ScfOptions,
    root: &Path,
) -> Result<Row<'a>> {
    let preopt = config.preopt.then(|| stage(start, &start.geometry, Surface::Gfnff, config, scf));
    let mut target = match &preopt {
        Some(preopt) => preopt
            .positions
            .as_ref()
            .map(|positions| stage(start, &start.geometry.with_bohr(positions), potential, config, scf)),
        None => Some(stage(start, &start.geometry, potential, config, scf)),
    };
    let stages = || preopt.iter().chain(&target);
    let steps = stages().map(|stage| stage.steps).sum();
    let evaluations = stages().map(|stage| stage.evaluations).sum();
    let hessian_evaluations = stages().map(|stage| stage.hessian_evaluations).sum();
    let potential_seconds = stages().map(|stage| stage.potential_seconds).sum();
    let setup_seconds = stages().map(|stage| stage.setup_seconds).sum();
    let gradient_seconds = stages().map(|stage| stage.gradient_seconds).sum();
    let hessian_seconds = stages().map(|stage| stage.hessian_seconds).sum();
    let optimizer_seconds = stages().map(|stage| stage.optimizer_seconds).sum();
    let wall_seconds = stages().map(|stage| stage.wall_seconds).sum();
    let reported = target.as_ref().or(preopt.as_ref()).context("no optimization stage")?;
    let scf_iterations_total = target.as_ref().and_then(|stage| stage.scf_iterations_total);
    let target_energy = target.as_ref().and_then(|stage| stage.energy);
    let max_gradient = target.as_ref().and_then(|stage| stage.max_gradient);
    let final_xyz = if let Some(positions) = target.as_ref().and_then(|stage| stage.positions.as_ref()) {
        let path = Path::new("geoms").join(potential.as_str()).join(&config.label).join(start.id.path());
        start.geometry.with_bohr(positions).write(&root.join(&path), &start.comment())?;
        Some(path)
    } else {
        None
    };
    let preopt_seconds = preopt.as_ref().map_or(0.0, |stage| stage.wall_seconds);
    let preopt_evaluations = preopt.as_ref().map_or(0, |stage| stage.evaluations);
    let trace = target
        .iter()
        .flat_map(|stage| &stage.trace)
        .map(|&Sample(n, seconds, energy, max_gradient, accuracy, gradient_norm)| {
            Sample(n + preopt_evaluations, seconds + preopt_seconds, energy, max_gradient, accuracy, gradient_norm)
        })
        .collect();
    let median_target_gradient = target.as_ref().and_then(|stage| stage.median_gradient_seconds);
    let median_gradient = if preopt.is_some() {
        let mut times: Vec<_> = stages().flat_map(|stage| stage.gradient_times.iter().copied()).collect();
        median_seconds(&mut times)
    } else {
        median_target_gradient
    };
    let preopt_gfn2_gradient_equivalents = if potential == Surface::Gfn2 {
        preopt.as_ref().zip(median_target_gradient).map(|(preopt, median)| preopt.gradient_seconds / median)
    } else {
        None
    };
    let total_gfn2_gradient_equivalents =
        preopt_gfn2_gradient_equivalents.map(|preopt| preopt + (evaluations - preopt_evaluations) as f64);
    let status = reported.status;
    let error = reported.error.clone();
    if let Some(stage) = &mut target {
        stage.trace.clear();
    }
    let row = Row {
        id: &start.id,
        potential,
        config: &config.label,
        atoms: start.geometry.numbers.len(),
        status,
        error,
        steps,
        evaluations,
        hessian_evaluations,
        potential_seconds,
        setup_seconds,
        gradient_seconds,
        gradient_evaluations: evaluations,
        hessian_seconds,
        hessian_gradient_equivalents: median_gradient.map(|seconds| hessian_seconds / seconds),
        optimizer_seconds,
        wall_seconds,
        scf_iterations_total,
        energy: target_energy,
        max_gradient,
        final_xyz,
        trace,
        preopt,
        target,
        preopt_gfn2_gradient_equivalents,
        total_gfn2_gradient_equivalents,
    };
    Ok(row)
}

fn status(status: &Status) -> (&'static str, Option<&str>) {
    match status {
        Status::Converged => ("converged", None),
        Status::MaxSteps => ("max_steps", None),
        Status::ScfFailure(error) => ("scf_failure", Some(error)),
        Status::Failed(error) => ("failed", Some(error)),
        Status::Cancelled => ("cancelled", None),
    }
}

fn stage(start: &Start, geometry: &Geometry, surface: Surface, config: &Config, scf: ScfOptions) -> Stage {
    let mut stage = Stage::new(surface);
    let beginning = Instant::now();
    let positions = geometry.bohr();
    let counter = Cell::new(0);
    let mut adapter = match BenchPotential::new(surface, start, geometry, config.hessian, scf, false, &counter) {
        Ok(adapter) => adapter,
        Err(error) => {
            stage.error = Some(format!("{error:#}"));
            stage.wall_seconds = beginning.elapsed().as_secs_f64();
            stage.setup_seconds = stage.wall_seconds;
            stage.potential_seconds = stage.setup_seconds;
            return stage;
        }
    };
    stage.setup_seconds = beginning.elapsed().as_secs_f64();
    let result = optimize(&mut adapter, &geometry.numbers, &positions, &config.options, |step| {
        stage.trace.push(Sample(
            counter.get(),
            beginning.elapsed().as_secs_f64(),
            step.energy,
            step.max_gradient,
            accuracy_name(step.accuracy),
            step.rms_gradient * ((3 * geometry.numbers.len()) as f64).sqrt(),
        ));
        ControlFlow::Continue(())
    });
    stage.wall_seconds = beginning.elapsed().as_secs_f64();
    let (status, error) = status(&result.status);
    stage.status = status;
    stage.error = error.map(str::to_owned);
    stage.method = Some(format!("{:?}", result.method));
    stage.steps = result.steps;
    stage.evaluations = result.evaluations;
    #[cfg(feature = "knobs")]
    {
        stage.hessian_evaluations = result.hessian_evaluations;
    }
    stage.gradient_seconds = adapter.calls.gradients.iter().sum();
    stage.potential_seconds = stage.setup_seconds + stage.gradient_seconds;
    stage.hessian_seconds = adapter.calls.hessian_seconds;
    stage.optimizer_seconds = (stage.wall_seconds - stage.potential_seconds - stage.hessian_seconds).max(0.0);
    stage.scf_iterations_total = (surface == Surface::Gfn2).then_some(adapter.calls.scf_iterations);
    if result.energy.is_finite() {
        stage.energy = Some(result.energy);
        stage.max_gradient = result.gradient.iter().map(|x| x.abs()).reduce(f64::max).filter(|value| value.is_finite());
        stage.positions = Some(result.positions);
    }
    stage.median_gradient_seconds = median_seconds(&mut adapter.calls.gradients);
    stage.gradient_times = adapter.calls.gradients;
    stage
}

fn median_seconds(gradients: &mut [f64]) -> Option<f64> {
    let size = gradients.len();
    if size == 0 {
        return None;
    }
    gradients.sort_unstable_by(f64::total_cmp);
    let median = (gradients[(size - 1) / 2] + gradients[size / 2]) * 0.5;
    (median > 0.0).then_some(median)
}

#[derive(clap::Args)]
pub struct ReferenceArgs {
    #[arg(long)]
    starts: PathBuf,
    #[arg(long)]
    potential: Surface,
    #[arg(long, required = true, num_args = 1..)]
    results: Vec<PathBuf>,
    #[arg(long)]
    out: PathBuf,
    #[arg(long, num_args = 1..)]
    xtb: Vec<PathBuf>,
}

#[derive(Deserialize)]
struct CandidateRow {
    #[serde(flatten)]
    id: Identity,
    potential: Surface,
    energy: Option<f64>,
    final_xyz: Option<PathBuf>,
}

#[derive(Serialize)]
struct ReferenceRow<'a> {
    #[serde(flatten)]
    id: &'a Identity,
    potential: Surface,
    energy: f64,
    xyz: PathBuf,
    status: &'static str,
    error: Option<String>,
}

pub fn reference(args: ReferenceArgs) -> Result<()> {
    let starts: BTreeMap<_, _> =
        starts(&args.starts, &Selection::default())?.into_iter().map(|start| (start.id.clone(), start)).collect();
    let mut candidates = BTreeMap::<Identity, Vec<PathBuf>>::new();
    for path in args.results.iter().chain(&args.xtb) {
        for (index, line) in
            BufReader::new(File::open(path).with_context(|| format!("opening {}", path.display()))?).lines().enumerate()
        {
            let row: CandidateRow =
                serde_json::from_str(&line?).with_context(|| format!("{}:{}", path.display(), index + 1))?;
            if row.potential == args.potential {
                ensure!(
                    starts.contains_key(&row.id),
                    "{} has no matching start in {}",
                    row.id.label(),
                    args.starts.display()
                );
                let paths = candidates.entry(row.id).or_default();
                if let (Some(_), Some(xyz)) = (row.energy, row.final_xyz) {
                    paths.push(parent(path).join(xyz));
                }
            }
        }
    }
    ensure!(!candidates.is_empty(), "no {} results to refine", args.potential.as_str());

    let root = parent(&args.out);
    fs::create_dir_all(root)?;
    let mut output = BufWriter::new(File::create(&args.out)?);
    let mut written = 0;
    for (id, paths) in candidates {
        let start = &starts[&id];
        if let Some(row) = refine(start, args.potential, paths, root)? {
            serde_json::to_writer(&mut output, &row)?;
            writeln!(output)?;
            written += 1;
        }
    }
    output.flush()?;
    ensure!(written > 0, "no usable reference geometries");
    println!("Wrote {written} references to {}", args.out.display());
    Ok(())
}

fn refine<'a>(
    start: &'a Start,
    surface: Surface,
    mut paths: Vec<PathBuf>,
    root: &Path,
) -> Result<Option<ReferenceRow<'a>>> {
    // Compare every candidate on the same surface, including official xtb output. In
    // particular, GFN-FF keeps the topology determined from the corresponding original start.
    let counter = Cell::new(0);
    let scf = ScfOptions::default();
    let mut adapter = match BenchPotential::new(
        surface,
        start,
        &start.geometry,
        Hessian::Model,
        scf.without_restart(),
        true,
        &counter,
    ) {
        Ok(adapter) => adapter,
        Err(error) => {
            eprintln!("{}: cannot build reference potential: {error:#}", start.id.label());
            return Ok(None);
        }
    };
    let mut best = None::<(f64, Vec<f64>)>;
    paths.sort_unstable();
    paths.dedup();
    for path in paths {
        let (geometry, _) = Geometry::read(&path)?;
        ensure!(geometry.numbers == start.geometry.numbers, "atom order differs in {}", path.display());
        let positions = geometry.bohr();
        let mut gradient = vec![0.0; positions.len()];
        match adapter.evaluate(&positions, &mut gradient, Accuracy::Tight) {
            Ok(energy) if energy.is_finite() => {
                if best.as_ref().is_none_or(|(current, _)| energy < *current) {
                    best = Some((energy, positions));
                }
            }
            result => eprintln!("{}: reference scoring failed: {result:?}", path.display()),
        }
    }
    let Some((mut energy, mut positions)) = best else {
        eprintln!("{}: no scored final geometry for reference", start.id.label());
        return Ok(None);
    };

    let geometry = start.geometry.with_bohr(&positions);
    let seed = positions.clone();
    let mut adapter = BenchPotential::new(surface, start, &geometry, Hessian::Model, scf, true, &counter)?;
    let options = Options { level: Level::Tight, max_steps: 3000, ..Options::default() };
    let result = optimize(&mut adapter, &geometry.numbers, &seed, &options, |step| {
        if step.energy.is_finite() && step.energy <= energy {
            energy = step.energy;
            positions.copy_from_slice(step.positions);
        }
        ControlFlow::Continue(())
    });
    let (status, error) = status(&result.status);
    if status != "converged" {
        eprintln!("{}: tight reference optimization ended with {:?}", start.id.label(), result.status);
    }
    let xyz = Path::new("references").join(surface.as_str()).join(start.id.path());
    start.geometry.with_bohr(&positions).write(&root.join(&xyz), &start.comment())?;
    Ok(Some(ReferenceRow { id: &start.id, potential: surface, energy, xyz, status, error: error.map(str::to_owned) }))
}
