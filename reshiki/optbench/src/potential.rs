// SPDX-License-Identifier: LGPL-3.0-or-later

use std::cell::Cell;
use std::time::Instant;

use anyhow::{Context, Result, ensure};
use axtb_driver::{Calculator, Error, Request, Task};
use clap::ValueEnum;
use cosmolkit::MolecularForceField;
use reshiki_optimize::{ANGSTROM_PER_BOHR, Accuracy, EvalError, KCAL_PER_HARTREE, Potential};

use crate::geometry::{Geometry, Start};

#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq, Serialize, ValueEnum)]
#[serde(rename_all = "lowercase")]
pub enum Surface {
    Gfn2,
    Gfnff,
    Mmff,
}
use serde::{Deserialize, Serialize};

impl Surface {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Gfn2 => "gfn2",
            Self::Gfnff => "gfnff",
            Self::Mmff => "mmff",
        }
    }

    fn calculator(self, start: &Start, geometry: &Geometry) -> Result<Calculator> {
        let mut request = Request::new(
            geometry.numbers.iter().map(|&z| usize::from(z)).collect(),
            geometry.bohr().as_chunks::<3>().0.to_vec(),
        );
        request.method = Some(self.as_str().into());
        request.charge = f64::from(start.charge);
        request.uhf = Some(i64::from(start.uhf));
        request.threads = Some(1);
        if self == Self::Gfnff {
            request.topology_positions = Some(start.geometry.bohr().as_chunks::<3>().0.to_vec());
        }
        Ok(Calculator::new(&request)?)
    }
}

#[derive(Clone, Copy, Default, PartialEq, Eq)]
pub enum Hessian {
    #[default]
    Model,
    Same,
    #[cfg(feature = "knobs")]
    Gfnff,
    #[cfg(feature = "knobs")]
    Auto,
}

#[derive(clap::Args, Clone, Copy)]
pub struct ScfOptions {
    #[arg(long)]
    no_restart: bool,
    #[arg(long, default_value_t = 30.0, value_parser = positive)]
    acc_coarse: f64,
    #[arg(long, default_value_t = 1.0, value_parser = positive)]
    acc_normal: f64,
    #[arg(long, default_value_t = 0.2, value_parser = positive)]
    acc_tight: f64,
}

fn positive(value: &str) -> Result<f64> {
    let value: f64 = value.parse()?;
    ensure!(value.is_finite() && value > 0.0, "accuracy must be positive and finite");
    Ok(value)
}

impl Default for ScfOptions {
    fn default() -> Self {
        Self { no_restart: false, acc_coarse: 30.0, acc_normal: 1.0, acc_tight: 0.2 }
    }
}

impl ScfOptions {
    pub fn without_restart(self) -> Self {
        Self { no_restart: true, ..self }
    }

    fn multiplier(self, accuracy: Accuracy) -> f64 {
        match accuracy_name(accuracy) {
            "normal" => self.acc_normal,
            "tight" => self.acc_tight,
            _ => self.acc_coarse,
        }
    }
}

pub fn accuracy_name(accuracy: Accuracy) -> &'static str {
    #[allow(unreachable_patterns)] // Older optimizer revisions have only Normal and Tight.
    match accuracy {
        Accuracy::Normal => "normal",
        Accuracy::Tight => "tight",
        _ => "coarse",
    }
}

enum Engine {
    Xtb(Box<Calculator>),
    Mmff(MolecularForceField),
}

#[derive(Default)]
pub struct Calls {
    pub gradients: Vec<f64>,
    pub hessian_seconds: f64,
    pub scf_iterations: usize,
}

pub struct BenchPotential<'a> {
    engine: Engine,
    surface: Surface,
    scf: ScfOptions,
    tight_scf: bool,
    counter: &'a Cell<usize>,
    pub calls: Calls,
    #[cfg(feature = "knobs")]
    hessian: Hessian,
    #[cfg(feature = "knobs")]
    gfnff: Option<Box<Calculator>>,
}

impl<'a> BenchPotential<'a> {
    pub fn new(
        surface: Surface,
        start: &Start,
        geometry: &Geometry,
        hessian: Hessian,
        scf: ScfOptions,
        tight_scf: bool,
        counter: &'a Cell<usize>,
    ) -> Result<Self> {
        ensure!(
            hessian != Hessian::Same || surface != Surface::Mmff,
            "COSMolKit MMFF94 does not provide an analytic Hessian"
        );
        let engine = match surface {
            Surface::Mmff => Engine::Mmff(start.mmff_field(geometry)?),
            _ => Engine::Xtb(Box::new(surface.calculator(start, geometry)?)),
        };
        Ok(Self {
            engine,
            surface,
            scf,
            tight_scf,
            counter,
            calls: Calls::default(),
            #[cfg(feature = "knobs")]
            hessian,
            #[cfg(feature = "knobs")]
            gfnff: if hessian == Hessian::Gfnff {
                Some(Box::new(Surface::Gfnff.calculator(start, geometry)?))
            } else {
                None
            },
        })
    }

    fn gradient(&mut self, positions: &[f64], gradient: &mut [f64], accuracy: Accuracy) -> Result<f64, EvalError> {
        let coordinates = positions.as_chunks::<3>().0;
        match &mut self.engine {
            Engine::Xtb(calculator) => {
                let task = Task { gradient: true, ..Task::default() };
                if self.surface == Surface::Gfn2 {
                    let accuracy = if self.tight_scf { Accuracy::Tight } else { accuracy };
                    calculator
                        .set_scf_accuracy(self.scf.multiplier(accuracy))
                        .map_err(|error| driver_error(error, &mut self.calls))?;
                    if self.scf.no_restart {
                        calculator.reset_wavefunction(false);
                    }
                }
                let response =
                    calculator.compute(coordinates, &task).map_err(|error| driver_error(error, &mut self.calls))?;
                count_scf(&response, &mut self.calls);
                let values = response.gradient.context("driver returned no gradient").map_err(failed)?;
                gradient.copy_from_slice(values.as_flattened());
                response.energy.context("driver returned no energy").map_err(failed)
            }
            Engine::Mmff(field) => {
                let angstrom: Vec<_> = coordinates.iter().map(|row| row.map(|x| x * ANGSTROM_PER_BOHR)).collect();
                field.set_positions_(&angstrom).map_err(failed)?;
                let evaluation = field.energy_gradient().map_err(failed)?;
                for (output, value) in gradient.iter_mut().zip(evaluation.gradient().as_flattened()) {
                    *output = value * ANGSTROM_PER_BOHR / KCAL_PER_HARTREE;
                }
                Ok(evaluation.energy() / KCAL_PER_HARTREE)
            }
        }
    }
}

impl Potential for BenchPotential<'_> {
    fn evaluate(&mut self, positions: &[f64], gradient: &mut [f64], accuracy: Accuracy) -> Result<f64, EvalError> {
        self.counter.set(self.counter.get() + 1);
        let start = Instant::now();
        let result = self.gradient(positions, gradient, accuracy);
        self.calls.gradients.push(start.elapsed().as_secs_f64());
        result
    }

    #[cfg(feature = "knobs")]
    fn hessian(&mut self, positions: &[f64], hessian: &mut [f64]) -> Option<Result<(), EvalError>> {
        if self.hessian == Hessian::Model {
            return None;
        }
        let calculator = match (&mut self.gfnff, &mut self.engine) {
            (Some(calculator), _) | (None, Engine::Xtb(calculator)) => calculator,
            (None, Engine::Mmff(_)) => return None,
        };
        let start = Instant::now();
        let result = (|| {
            if self.scf.no_restart && calculator.method() == axtb_driver::Method::Gfn2 {
                calculator.reset_wavefunction(false);
            }
            let request = Request {
                hessian: true,
                hessian_method: Some("analytic".into()),
                vibrations: Some(false),
                ..Request::default()
            };
            let task = Task::from_request(&request).map_err(failed)?;
            let response = calculator
                .compute(positions.as_chunks::<3>().0, &task)
                .map_err(|error| driver_error(error, &mut self.calls))?;
            count_scf(&response, &mut self.calls);
            let matrix = response.hessian.context("driver returned no Hessian").map_err(failed)?.hessian;
            ensure_hessian_size(&matrix, positions.len(), hessian.len())?;
            hessian.copy_from_slice(&matrix.data);
            Ok(())
        })();
        self.calls.hessian_seconds += start.elapsed().as_secs_f64();
        Some(result)
    }
}

#[cfg(feature = "knobs")]
fn ensure_hessian_size(matrix: &axtb_driver::Matrix, dimension: usize, size: usize) -> Result<(), EvalError> {
    if matrix.rows == dimension && matrix.cols == dimension && matrix.data.len() == size {
        Ok(())
    } else {
        Err(failed("driver returned a Hessian with the wrong size"))
    }
}

fn count_scf(response: &axtb_driver::Response, calls: &mut Calls) {
    if let Some(scf) = &response.scf {
        calls.scf_iterations += scf
            .scf_attempts
            .as_ref()
            .filter(|attempts| !attempts.is_empty())
            .map_or(scf.iterations, |attempts| attempts.iter().filter_map(|a| a.iterations).sum());
    }
}

fn driver_error(error: Error, calls: &mut Calls) -> EvalError {
    if let Error::ScfNotConverged { iterations, .. } = error {
        calls.scf_iterations += iterations;
        EvalError::Scf(error.to_string())
    } else {
        failed(error)
    }
}

fn failed(error: impl std::fmt::Display) -> EvalError {
    EvalError::Failed(error.to_string())
}
