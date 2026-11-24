// SPDX-License-Identifier: LGPL-3.0-or-later
//! Names of methods, conventions and settings as every front end accepts them.
//!
//! One alias rule for all of them: case and everything but ASCII letters and digits are
//! ignored, so `"GFN2-xTB"`, `"gfn2_xtb"` and `"gfn2xtb"` are one name, as are `"xtb-6.7.1"`
//! and `"XTB 6.7.1"`. The responses report the canonical spelling.

use crate::error::{Error, Result};
use axtb::vibrations::Conventions as VibConventions;

/// The key of a name under the alias rule: its ASCII letters and digits, lower case.
pub fn key(name: &str) -> String {
    name.chars().filter(char::is_ascii_alphanumeric).map(|c| c.to_ascii_lowercase()).collect()
}

/// The methods of the front ends.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Method {
    /// GFN2-xTB (`gfn2`).
    Gfn2,
    /// GFN-FF (`gfnff`).
    Gfnff,
    /// g-xTB v1.1 (`gxtb-v1.1`), the public preprint version.
    GxtbV1_1,
    /// g-xTB 2.0.1 (`gxtb-v2.0.1`).
    GxtbV2_0_1,
}

impl Method {
    /// The method a name selects: `gfn2` (also `GFN2-xTB`), `gfnff`, `gxtb-v1.1` (also
    /// `gxtb-1.1`, `g-xTB 1.1`) or `gxtb-v2.0.1` (also `gxtb-2.0.1`, `gxtb-v2.0`, `g-xTB 2.0`).
    pub fn parse(name: &str) -> Result<Self> {
        match key(name).as_str() {
            "gfn2" | "gfn2xtb" => Ok(Method::Gfn2),
            "gfnff" => Ok(Method::Gfnff),
            "gxtbv11" | "gxtb11" => Ok(Method::GxtbV1_1),
            "gxtbv201" | "gxtb201" | "gxtbv20" | "gxtb20" => Ok(Method::GxtbV2_0_1),
            _ => Err(Error::Input(format!(
                "unknown method {name:?}; use \"gfn2\", \"gfnff\", \"gxtb-v1.1\" or \"gxtb-v2.0.1\""
            ))),
        }
    }

    /// The canonical name: `gfn2`, `gfnff`, `gxtb-v1.1` or `gxtb-v2.0.1`.
    pub fn name(self) -> &'static str {
        match self {
            Method::Gfn2 => "gfn2",
            Method::Gfnff => "gfnff",
            Method::GxtbV1_1 => "gxtb-v1.1",
            Method::GxtbV2_0_1 => "gxtb-v2.0.1",
        }
    }

    /// Whether the method has an SCF (every method but GFN-FF).
    pub fn is_scf(self) -> bool {
        self != Method::Gfnff
    }

    /// Whether the method is a g-xTB version.
    pub fn is_gxtb(self) -> bool {
        matches!(self, Method::GxtbV1_1 | Method::GxtbV2_0_1)
    }

    /// Check that this build includes the method. A disabled method is an input error.
    pub fn check_available(self) -> Result<()> {
        if !cfg!(feature = "gxtb") && self.is_gxtb() {
            return Err(Error::gxtb_unavailable());
        }
        Ok(())
    }

    /// The g-xTB version of a g-xTB method.
    #[cfg(feature = "gxtb")]
    pub(crate) fn gxtb_version(self) -> Option<axtb_gxtb::Version> {
        match self {
            Method::GxtbV1_1 => Some(axtb_gxtb::Version::V1_1),
            Method::GxtbV2_0_1 => Some(axtb_gxtb::Version::V2_0_1),
            _ => None,
        }
    }
}

impl std::fmt::Display for Method {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.name())
    }
}

/// The canonical name of a method name (see [`Method::parse`]); the codec and Python use it.
pub fn canonical_method(name: &str) -> Result<&'static str> {
    Method::parse(name).map(Method::name)
}

/// GFN2-xTB conventions: `axtb` (also `default`), `xtb` (also `xtb-6.8.0`), `xtb-6.7.1` or
/// `tblite`.
pub fn gfn2_conventions(name: &str) -> Result<axtb::Conventions> {
    match key(name).as_str() {
        "axtb" | "default" => Ok(axtb::Conventions::Axtb),
        "xtb" | "xtb680" => Ok(axtb::Conventions::Xtb),
        "xtb671" => Ok(axtb::Conventions::Xtb671),
        "tblite" => Ok(axtb::Conventions::Tblite),
        _ => Err(Error::Input(format!(
            "unknown gfn2 conventions {name:?}; use \"axtb\", \"xtb\" (xtb 6.8.0), \"xtb-6.7.1\" or \"tblite\""
        ))),
    }
}

/// g-xTB conventions: `axtb` (also `default`) or `gxtb` (the program's).
#[cfg(feature = "gxtb")]
pub(crate) fn gxtb_conventions(name: &str) -> Result<axtb_gxtb::Conventions> {
    match key(name).as_str() {
        "axtb" | "default" => Ok(axtb_gxtb::Conventions::Axtb),
        "gxtb" | "gxtbprogram" => Ok(axtb_gxtb::Conventions::Gxtb),
        _ => Err(Error::Input(format!(
            "unknown g-xTB conventions {name:?}; use \"axtb\" or \"gxtb\" (the g-xTB program's)"
        ))),
    }
}

/// GFN-FF conventions: `axtb` (also `default`; xtb 6.8.0's force field), `xtb` (also
/// `xtb-6.8.0`) or `xtb-6.7.1`.
pub(crate) fn gfnff_conventions(name: &str) -> Result<axtb_gfnff::Conventions> {
    use axtb_gfnff::Conventions as Ff;
    match key(name).as_str() {
        "axtb" | "default" => Ok(Ff::Axtb),
        "xtb" | "xtb680" => Ok(Ff::Xtb680),
        "xtb671" => Ok(Ff::Xtb671),
        _ => Err(Error::Input(format!(
            "unknown gfnff conventions {name:?}; use \"axtb\", \"xtb\" (xtb 6.8.0) or \"xtb-6.7.1\""
        ))),
    }
}

/// Conventions of the vibrational analysis: `axtb` (also `default`), `xtb` (also
/// `xtb-6.8.0`), `xtb-6.7.1` (also `gxtb`: the g-xTB 2.0.1 program is a 6.7.1-based xtb).
pub fn analysis_conventions(name: &str) -> Result<VibConventions> {
    match key(name).as_str() {
        "axtb" | "default" => Ok(VibConventions::Axtb),
        "xtb" | "xtb680" => Ok(VibConventions::Xtb),
        "xtb671" | "gxtb" => Ok(VibConventions::Xtb671),
        _ => Err(Error::Input(format!(
            "unknown vibrational-analysis conventions {name:?}; use \"axtb\", \"xtb\" (xtb 6.8.0) or \"xtb-6.7.1\""
        ))),
    }
}

/// Hessian method: `analytic`, `response` or `numerical` (also `numeric`, `fd`).
pub fn hessian_method(name: &str) -> Result<axtb::HessianMethod> {
    match key(name).as_str() {
        "analytic" | "analytical" => Ok(axtb::HessianMethod::Analytic),
        "response" => Ok(axtb::HessianMethod::Response),
        "numerical" | "numeric" | "fd" => Ok(axtb::HessianMethod::Numerical),
        _ => Err(Error::Input(format!(
            "unknown Hessian method {name:?}; use \"analytic\", \"response\" or \"numerical\""
        ))),
    }
}

/// The canonical name of a Hessian method.
pub fn hessian_method_label(m: axtb::HessianMethod) -> &'static str {
    match m {
        axtb::HessianMethod::Analytic => "analytic",
        axtb::HessianMethod::Response => "response",
        axtb::HessianMethod::Numerical => "numerical",
    }
}

/// GFN2-xTB SCF mixer: `auto`, `broyden` or `diis` (also `pulay`).
pub fn mixer(name: &str) -> Result<axtb::MixerKind> {
    match key(name).as_str() {
        "auto" => Ok(axtb::MixerKind::Auto),
        "broyden" => Ok(axtb::MixerKind::Broyden),
        "diis" | "pulay" => Ok(axtb::MixerKind::Diis),
        _ => Err(Error::Input(format!("unknown SCF mixer {name:?}; use \"auto\", \"broyden\" or \"diis\""))),
    }
}

/// Reference state of the solvation free energy: `gsolv`, `bar1M` or `reference`.
pub fn solvation_state(name: &str) -> Result<axtb::ReferenceState> {
    match key(name).as_str() {
        "gsolv" => Ok(axtb::ReferenceState::Gsolv),
        "bar1m" => Ok(axtb::ReferenceState::Bar1M),
        "reference" => Ok(axtb::ReferenceState::Reference),
        _ => Err(Error::Input(format!(
            "unknown solvation reference state {name:?}; use \"gsolv\", \"bar1M\" or \"reference\""
        ))),
    }
}

/// Functions of a Molden file: `spherical` (also `pure`, `5d`) or `cartesian` (also `6d`).
pub fn molden_functions(name: &str) -> Result<axtb::MoldenFunctions> {
    match key(name).as_str() {
        "spherical" | "pure" | "5d" => Ok(axtb::MoldenFunctions::Spherical),
        "cartesian" | "6d" => Ok(axtb::MoldenFunctions::Cartesian),
        _ => Err(Error::Input(format!("unknown Molden functions {name:?}; use \"spherical\" or \"cartesian\""))),
    }
}

/// Solvation correction of the g-xTB composite: `gfn2` (also `GFN2-xTB`).
#[cfg(feature = "gxtb")]
pub(crate) fn correction(name: &str) -> Result<axtb_gxtb::Correction> {
    match key(name).as_str() {
        "gfn2" | "gfn2xtb" => Ok(axtb_gxtb::Correction::Gfn2),
        _ => Err(Error::Input(format!(
            "unknown solvation correction {name:?}; the only one is \"gfn2\" (GFN2-xTB with ALPB/GBSA)"
        ))),
    }
}

/// The restart policy of a calculator's SCFs (`restart`): `previous` (also `restart`) or
/// `fresh` (also `off`).
pub fn restart(name: &str) -> Result<axtb::method::Restart> {
    match key(name).as_str() {
        "previous" | "restart" => Ok(axtb::method::Restart::Previous),
        "fresh" | "off" => Ok(axtb::method::Restart::Fresh),
        _ => Err(Error::Input(format!("unknown restart {name:?}; use \"previous\" or \"fresh\""))),
    }
}

/// The canonical name of `name` of a kind (for the Python binding): `method` (see
/// [`Method::parse`]); `conventions` of `method` (required; e.g. `xtb-6.8.0` -> `xtb`); `solvent`
/// (ALPB's table name, else GBSA's; e.g. `H2O` -> `water`).
pub fn canonical(kind: &str, name: &str, method: Option<&str>) -> Result<&'static str> {
    match key(kind).as_str() {
        "method" => canonical_method(name),
        "conventions" => {
            let m = method.ok_or_else(|| Error::Input("canonical conventions need the method".into()))?;
            Ok(match Method::parse(m)? {
                Method::Gfn2 => gfn2_conventions(name)?.label(),
                Method::Gfnff => gfnff_conventions(name)?.label(),
                #[cfg(feature = "gxtb")]
                Method::GxtbV1_1 | Method::GxtbV2_0_1 => gxtb_conventions(name)?.label(),
                #[cfg(not(feature = "gxtb"))]
                Method::GxtbV1_1 | Method::GxtbV2_0_1 => return Err(Error::gxtb_unavailable()),
            })
        }
        "solvent" => [axtb_solvation::Model::Alpb, axtb_solvation::Model::Gbsa]
            .into_iter()
            .find_map(|m| axtb_solvation::canonical_solvent(name, m))
            .ok_or_else(|| Error::Input(format!("unknown solvent {name:?}"))),
        _ => {
            Err(Error::Input(format!("unknown kind of name {kind:?}; use \"method\", \"conventions\" or \"solvent\"")))
        }
    }
}

/// The electrons a neutral atom of atomic number `z` brings to the orbitals of `method`: the sum
/// of its reference shell occupations (the shells of its basis), so Z less them are its core
/// electrons (`Electronic::core_electrons`); GFN-FF, which has no orbitals, all Z. `None` where the
/// method is unavailable or its tables end (GFN2-xTB above Rn, g-xTB above Lr; [`supported_elements`] says which
/// elements a method has parameters for). For the Python binding's check of the spin parity
/// before any calculation runs.
pub fn valence_electrons(method: Method, z: usize) -> Option<f64> {
    if z == 0 {
        return None;
    }
    let i = z - 1;
    match method {
        Method::Gfnff => Some(z as f64),
        Method::Gfn2 => {
            use axtb::param::gfn2_data::{ANG_SHELL, NSHELL, REFERENCE_OCC};
            let n = *NSHELL.get(i)?;
            Some((0..n).map(|ish| REFERENCE_OCC[i][ANG_SHELL[i][ish]]).sum())
        }
        #[cfg(feature = "gxtb")]
        Method::GxtbV2_0_1 => axtb_gxtb::data::v2_0_1::REF_OCC.get(i).map(|occ| occ.iter().sum()),
        #[cfg(feature = "gxtb")]
        Method::GxtbV1_1 => axtb_gxtb::data::runtime::REF_SHELL_OCC.get(i).map(|occ| occ.iter().sum()),
        #[cfg(not(feature = "gxtb"))]
        Method::GxtbV1_1 | Method::GxtbV2_0_1 => None,
    }
}

/// One shell of a method's basis for an element ([`shell_layout`]).
#[derive(Clone, Debug, PartialEq)]
pub struct ShellLayout {
    /// The atom (index into the atomic numbers).
    pub atom: usize,
    /// Angular momentum.
    pub l: usize,
    /// The reference occupation the shell's Mulliken charge is relative to.
    pub reference_occupation: f64,
    /// The primitive exponents where they do not depend on the geometry (g-xTB's q-vSZP); empty
    /// for GFN2-xTB.
    pub exponents: Vec<f64>,
}

/// The shells of `method`'s basis for the atomic numbers `numbers`, atom by atom in the basis'
/// order, from the parameter tables (no calculation): for the Python binding's check of a given
/// wavefunction or given moments against the method. GFN-FF has no shells. Returns `None` for
/// unavailable methods and atomic numbers beyond the tables.
pub fn shell_layout(method: Method, numbers: &[usize]) -> Option<Vec<ShellLayout>> {
    method.check_available().ok()?;
    let mut out = Vec::new();
    for (atom, &z) in numbers.iter().enumerate() {
        let i = z.checked_sub(1)?;
        match method {
            Method::Gfnff => return Some(Vec::new()),
            Method::Gfn2 => {
                use axtb::param::gfn2_data::{ANG_SHELL, NSHELL, REFERENCE_OCC};
                for ish in 0..*NSHELL.get(i)? {
                    let l = ANG_SHELL[i][ish];
                    out.push(ShellLayout { atom, l, reference_occupation: REFERENCE_OCC[i][l], exponents: Vec::new() });
                }
            }
            #[cfg(feature = "gxtb")]
            Method::GxtbV1_1 | Method::GxtbV2_0_1 => {
                let element = axtb_gxtb::data::qvszp::ELEMENTS.get(z)?.as_ref()?;
                let reference: &[f64; 4] = if method == Method::GxtbV2_0_1 {
                    axtb_gxtb::data::v2_0_1::REF_OCC.get(i)?
                } else {
                    axtb_gxtb::data::runtime::REF_SHELL_OCC.get(i)?
                };
                for (ish, sh) in element.shells.iter().enumerate() {
                    out.push(ShellLayout {
                        atom,
                        l: sh.l,
                        reference_occupation: reference[ish],
                        exponents: sh.alpha.to_vec(),
                    });
                }
            }
            #[cfg(not(feature = "gxtb"))]
            Method::GxtbV1_1 | Method::GxtbV2_0_1 => return None,
        }
    }
    Some(out)
}

/// The atomic numbers `method` has parameters for (every element a one-atom calculation of the
/// method can be set up for; computed once per process). A disabled method has no supported elements.
pub fn supported_elements(method: Method) -> Vec<usize> {
    if method.check_available().is_err() {
        return Vec::new();
    }
    use std::sync::OnceLock;
    static TABLE: [OnceLock<Vec<usize>>; 4] = [OnceLock::new(), OnceLock::new(), OnceLock::new(), OnceLock::new()];
    let i = match method {
        Method::Gfn2 => 0,
        Method::Gfnff => 1,
        Method::GxtbV1_1 => 2,
        Method::GxtbV2_0_1 => 3,
    };
    TABLE[i]
        .get_or_init(|| {
            (1..=118)
                .filter(|&z| {
                    let req = crate::Request {
                        method: Some(method.name().into()),
                        ..crate::Request::new(vec![z], vec![[0.0; 3]])
                    };
                    !matches!(crate::Calculator::new(&req), Err(Error::UnsupportedElement { .. }))
                })
                .collect()
        })
        .clone()
}

/// Where a calculator's SCFs start ([`guess`]).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Guess {
    /// The previous geometry's converged state (`previous`, also `restart`).
    Previous,
    /// The method's own start at every geometry (`fresh`).
    Fresh,
    /// Every geometry from EEQ charges (`eeq`; GFN2-xTB's own start).
    Eeq,
    /// Every geometry from neutral atoms (`sad`; GFN2-xTB).
    Sad,
    /// Every geometry from the program's EEQ-BC start (`eeq-bc`; g-xTB's own start).
    EeqBc,
}

/// Where a calculator's SCFs start: `previous` (also `restart`), the previous geometry's
/// converged state; `fresh`, the method's own start; or a named start, `eeq` or `sad`
/// (GFN2-xTB) or `eeq-bc` (g-xTB), which the calculator checks against its method.
pub fn guess(name: &str) -> Result<Guess> {
    match key(name).as_str() {
        "previous" | "restart" => Ok(Guess::Previous),
        "fresh" => Ok(Guess::Fresh),
        "eeq" => Ok(Guess::Eeq),
        "sad" => Ok(Guess::Sad),
        "eeqbc" => Ok(Guess::EeqBc),
        _ => Err(Error::Input(format!(
            "unknown guess {name:?}; use \"previous\", \"fresh\" (the method's own start), \"eeq\" or \"sad\" \
             (GFN2-xTB) or \"eeq-bc\" (g-xTB)"
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn one_alias_rule() {
        for (name, method) in [
            ("gfn2", "gfn2"),
            ("GFN2-xTB", "gfn2"),
            ("GFN2 xTB", "gfn2"),
            ("gfn-ff", "gfnff"),
            ("g-xTB 1.1", "gxtb-v1.1"),
            ("gxtb_v1.1", "gxtb-v1.1"),
            ("gxtb-2.0.1", "gxtb-v2.0.1"),
            ("g-xTB 2.0", "gxtb-v2.0.1"),
        ] {
            assert_eq!(canonical_method(name).unwrap(), method, "{name}");
        }
        assert!(matches!(canonical_method("gfn1"), Err(Error::Input(_))));
        assert_eq!(gfn2_conventions("XTB 6.7.1").unwrap(), axtb::Conventions::Xtb671);
        assert_eq!(gfn2_conventions("default").unwrap(), axtb::Conventions::Axtb);
        assert_eq!(gfn2_conventions("xtb-6.8.0").unwrap(), axtb::Conventions::Xtb);
        assert_eq!(gfnff_conventions("xtb_6.7.1").unwrap().label(), "xtb-6.7.1");
        assert_eq!(analysis_conventions("gxtb").unwrap(), VibConventions::Xtb671);
        assert_eq!(solvation_state("BAR 1M").unwrap(), axtb::ReferenceState::Bar1M);
        assert_eq!(guess("eeq-bc").unwrap(), Guess::EeqBc);
        assert_eq!(guess("Restart").unwrap(), Guess::Previous);
        assert_eq!(guess("SAD").unwrap(), Guess::Sad);
        assert!(guess("huckel").is_err());
        assert!(hessian_method("finite").is_err());
    }
}
