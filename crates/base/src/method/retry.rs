// SPDX-License-Identifier: LGPL-3.0-or-later
//! The one retry policy: a ladder of tries climbed until one succeeds.
//!
//! A caller's SCF retries ([`Retry`], e.g. axtb-driver's `scf_retries`) become fresh rungs after
//! the first try ([`Retry::rungs`]); every try that ran can be logged ([`Try`]).

use super::{Method, ScfSettings};
use crate::{Error, Result};

/// One try of a calculation: what changes against the method's settings.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Rung {
    /// Forget the wavefunction before this try: start from the method's own start.
    pub fresh: bool,
    /// Electronic temperature (K) of this try; `None`: the method's setting.
    pub etemp: Option<f64>,
    /// Iteration limit of this try; `None`: the method's setting.
    pub max_iter: Option<usize>,
    /// How reports name this try (`None`: the first, unchanged one).
    pub label: Option<String>,
}

impl Rung {
    /// The SCF settings of this try: `base` with the rung's electronic temperature and iteration
    /// limit where it has them.
    pub fn apply(&self, base: &ScfSettings) -> ScfSettings {
        let mut s = base.clone();
        s.etemp = self.etemp.or(base.etemp);
        s.max_iter = self.max_iter.or(base.max_iter);
        s
    }
}

/// One SCF retry a caller asks for: a fresh SCF with another electronic temperature (K) or
/// iteration limit; a missing value keeps the calculation's.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Retry {
    pub etemp: Option<f64>,
    pub max_iter: Option<usize>,
}

impl Retry {
    /// Refuse retries that cannot run or change nothing: none at all, a rung without a value, a
    /// temperature that is negative, not finite or below `floor` (K: an explicit electronic
    /// temperature of the calculation, which no rung may go below), an iteration limit of 0.
    /// The messages name the rung as `scf_retries[i]`.
    pub fn check(retries: &[Retry], floor: Option<f64>) -> Result<()> {
        if retries.is_empty() {
            return Err(Error::Input("scf_retries: at least one rung (or leave scf_retries out)".into()));
        }
        for (i, r) in retries.iter().enumerate() {
            if r.etemp.is_none() && r.max_iter.is_none() {
                return Err(Error::Input(format!(
                    "scf_retries[{i}]: a rung changes etemp or max_iter (one without either repeats the SCF that \
                     failed)"
                )));
            }
            if let Some(t) = r.etemp {
                if !(t >= 0.0 && t.is_finite()) {
                    return Err(Error::Input(format!("scf_retries[{i}].etemp: must be >= 0 K and finite, got {t}")));
                }
                if let Some(f) = floor.filter(|&f| t < f) {
                    return Err(Error::Input(format!(
                        "scf_retries[{i}].etemp: {t} K is below the calculation's electronic temperature ({f} K), \
                         which the retries start from"
                    )));
                }
            }
            if r.max_iter == Some(0) {
                return Err(Error::Input(format!("scf_retries[{i}].max_iter: must be at least 1")));
            }
        }
        Ok(())
    }

    /// The rungs of `retries` after a first try at `base` (electronic temperature, iteration
    /// limit): each from the method's own start, with its values (a missing one: the base's) and
    /// the label `retry i of n: etemp T K, max_iter N`.
    pub fn rungs(base: (f64, usize), retries: &[Retry]) -> Vec<Rung> {
        let n = retries.len();
        retries
            .iter()
            .enumerate()
            .map(|(i, r)| {
                let (t, m) = (r.etemp.unwrap_or(base.0), r.max_iter.unwrap_or(base.1));
                Rung {
                    fresh: true,
                    etemp: Some(t),
                    max_iter: Some(m),
                    label: Some(format!("retry {} of {n}: etemp {t} K, max_iter {m}", i + 1)),
                }
            })
            .collect()
    }
}

/// One try that ran on a retry ladder: the part, the rung and its settings, how it ended.
#[derive(Clone, Debug, PartialEq)]
pub struct Try {
    /// The part (a composite's key); `None`: the method itself.
    pub part: Option<&'static str>,
    /// The rung: 0 the first try, then the ladder's rungs in order.
    pub rung: usize,
    /// Electronic temperature (K) and iteration limit of the try.
    pub etemp: f64,
    pub max_iter: usize,
    /// Where its SCF started and its iterations, where known (a single point's; an SCF that did
    /// not converge: its iterations).
    pub start: Option<&'static str>,
    pub iterations: Option<usize>,
    /// The mixer that converged it (GFN2-xTB).
    pub mixer: Option<&'static str>,
    /// `None`: it succeeded; else its error.
    pub error: Option<Error>,
}

/// Which errors climb to the next rung; every other error ends the climb.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RetryOn {
    /// An SCF that did not converge ([`Error::ScfNotConverged`]).
    NotConverged,
    /// An SCF that failed: not converged or broken down ([`Error::Numerical`], e.g. a failed
    /// mixer or a lost Fermi level).
    Failed,
}

impl RetryOn {
    /// Whether `e` climbs.
    pub fn retries(self, e: &Error) -> bool {
        match self {
            RetryOn::NotConverged => matches!(e, Error::ScfNotConverged { .. }),
            RetryOn::Failed => matches!(e, Error::ScfNotConverged { .. } | Error::Numerical(_)),
        }
    }
}

/// Tries in order ([`RetryLadder::climb`]).
#[derive(Clone, Debug, PartialEq)]
pub struct RetryLadder {
    pub rungs: Vec<Rung>,
    pub on: RetryOn,
}

impl RetryLadder {
    /// One try, as the method stands.
    pub fn once() -> Self {
        Self { rungs: vec![Rung::default()], on: RetryOn::NotConverged }
    }

    /// The SCF as the method stands, and when `warm` (it started from a previous wavefunction)
    /// and did not converge, once more from the method's own start.
    pub fn warm_then_fresh(warm: bool) -> Self {
        let mut rungs = vec![Rung::default()];
        if warm {
            rungs.push(Rung { fresh: true, label: Some("fresh start".into()), ..Rung::default() });
        }
        Self { rungs, on: RetryOn::NotConverged }
    }

    /// Run `attempt` with each rung in order until one succeeds; an error the ladder does not
    /// retry ends the climb with it, and so does the last rung's error. Returns the result and
    /// the index of the rung that gave it.
    pub fn climb<T>(&self, mut attempt: impl FnMut(&Rung) -> Result<T>) -> Result<(T, usize)> {
        let mut last = None;
        for (i, rung) in self.rungs.iter().enumerate() {
            match attempt(rung) {
                Ok(v) => return Ok((v, i)),
                Err(e) if self.on.retries(&e) => last = Some(e),
                Err(e) => return Err(e),
            }
        }
        Err(last.unwrap_or_else(|| Error::Input("a retry ladder without rungs".into())))
    }

    /// [`RetryLadder::climb`] with `methods` tried together at each rung, from the rung `from`
    /// on: at each rung every method in order at the rung's settings over its own
    /// ([`Rung::apply`]; a method without SCF settings as it is), its wavefunction forgotten
    /// first on a fresh rung, and `f(k, method, i, rung)` tries method `k` at rung `i`; the first
    /// that fails ends the rung's try. Afterwards every method has its own settings again.
    /// Returns each method's result and the index of the rung that gave them.
    pub fn climb_methods<T>(
        &self,
        methods: &mut [&mut dyn Method],
        from: usize,
        mut f: impl FnMut(usize, &mut dyn Method, usize, &Rung) -> Result<T>,
    ) -> Result<(Vec<T>, usize)> {
        let base: Vec<_> = methods.iter().map(|m| m.scf_settings().cloned()).collect();
        let ladder = RetryLadder { rungs: self.rungs.get(from..).unwrap_or_default().to_vec(), on: self.on };
        let mut index = from;
        let r = ladder.climb(|rung| {
            let i = index;
            index += 1;
            let mut results = Vec::with_capacity(methods.len());
            for (k, (m, b)) in methods.iter_mut().zip(&base).enumerate() {
                if let Some(b) = b {
                    m.set_scf_settings(rung.apply(b))?;
                }
                if rung.fresh {
                    m.reset_wavefunction();
                }
                results.push(f(k, &mut **m, i, rung)?);
            }
            Ok(results)
        });
        for (m, b) in methods.iter_mut().zip(base) {
            if let Some(b) = b {
                m.set_scf_settings(b)?;
            }
        }
        r.map(|(v, i)| (v, i + from))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn nc(i: usize) -> Error {
        Error::ScfNotConverged { iterations: i, energy_change: 0.0, residual: 0.0 }
    }

    #[test]
    fn climbs_until_one_succeeds_and_stops_on_other_errors() {
        let ladder = RetryLadder::warm_then_fresh(true);
        let mut tried = Vec::new();
        let r = ladder.climb(|rung| {
            tried.push(rung.fresh);
            if rung.fresh {
                Ok(7)
            } else {
                Err(nc(1))
            }
        });
        assert_eq!(r, Ok((7, 1)));
        assert_eq!(tried, [false, true]);
        // the last rung's error; a numerical failure ends a NotConverged ladder at once
        assert_eq!(ladder.climb(|r| Err::<(), _>(nc(r.fresh as usize))), Err(nc(1)));
        let mut n = 0;
        let e = ladder.climb(|_| {
            n += 1;
            Err::<(), _>(Error::Numerical("x".into()))
        });
        assert_eq!((e, n), (Err(Error::Numerical("x".into())), 1));
        let failed = RetryLadder { on: RetryOn::Failed, ..ladder };
        assert_eq!(failed.climb(|r| if r.fresh { Ok(()) } else { Err(Error::Numerical("x".into())) }), Ok(((), 1)));
        assert_eq!(RetryLadder::warm_then_fresh(false).rungs.len(), 1);
    }

    #[test]
    fn retries_become_fresh_labelled_rungs_and_are_checked() {
        let given = [Retry { etemp: Some(1000.0), max_iter: Some(1000) }, Retry { etemp: None, max_iter: Some(500) }];
        let rungs = Retry::rungs((300.0, 250), &given);
        assert_eq!(
            rungs[0],
            Rung {
                fresh: true,
                etemp: Some(1000.0),
                max_iter: Some(1000),
                label: Some("retry 1 of 2: etemp 1000 K, max_iter 1000".into())
            }
        );
        assert_eq!(
            (rungs[1].etemp, rungs[1].label.as_deref()),
            (Some(300.0), Some("retry 2 of 2: etemp 300 K, max_iter 500"))
        );
        // a rung's settings over the method's: a missing value keeps the method's
        let base = ScfSettings { accuracy: 0.5, etemp: None, max_iter: Some(40), ..ScfSettings::default() };
        let r = Rung { max_iter: Some(9), ..Rung::default() }.apply(&base);
        assert_eq!((r.accuracy, r.etemp, r.max_iter), (0.5, None, Some(9)));
        assert_eq!(Rung::default().apply(&base), base);
        // the checks
        assert!(Retry::check(&given, None).is_ok() && Retry::check(&given, Some(1000.0)).is_ok());
        let floor = Retry::check(&given, Some(1500.0)).unwrap_err().to_string();
        assert!(floor.contains("scf_retries[0].etemp") && floor.contains("below"), "{floor}");
        for bad in [
            vec![],
            vec![Retry::default()],
            vec![Retry { etemp: Some(-1.0), max_iter: None }],
            vec![Retry { etemp: Some(f64::NAN), max_iter: None }],
            vec![Retry { etemp: None, max_iter: Some(0) }],
        ] {
            assert!(matches!(Retry::check(&bad, None), Err(Error::Input(_))), "{bad:?}");
        }
    }
}
