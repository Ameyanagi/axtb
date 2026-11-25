// SPDX-License-Identifier: LGPL-3.0-or-later

use std::collections::HashSet;
use std::fmt::Write as _;
use std::fs;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail, ensure};
use cosmolkit::{Element, EmbedParams, ForceFieldMinimizeParams, MmffForceFieldParams, MolecularForceField, Molecule};
use reshiki_optimize::ANGSTROM_PER_BOHR;
use serde::{Deserialize, Serialize};

#[derive(Clone)]
pub struct Geometry {
    pub numbers: Vec<u8>,
    pub angstrom: Vec<[f64; 3]>,
}

impl Geometry {
    pub fn read(path: &Path) -> Result<(Self, String)> {
        let read = || {
            let file = fs::read_to_string(path)?;
            let mut lines = file.lines();
            let atoms: usize = lines.next().context("missing atom count")?.trim().parse()?;
            ensure!(atoms > 0, "empty XYZ");
            let comment = lines.next().context("missing XYZ comment")?.to_owned();
            let mut numbers = Vec::new();
            let mut angstrom = Vec::new();
            for line in lines.by_ref().take(atoms) {
                let mut fields = line.split_whitespace();
                let symbol = fields.next().context("missing element")?;
                let element = Element::from_symbol(symbol)
                    .or_else(|| symbol.parse().ok().and_then(Element::from_atomic_number))
                    .or_else(|| {
                        let mut normalized = symbol.to_ascii_lowercase().into_bytes();
                        normalized.first_mut()?.make_ascii_uppercase();
                        Element::from_symbol(std::str::from_utf8(&normalized).ok()?)
                    })
                    .with_context(|| format!("invalid element {symbol:?}"))?;
                ensure!(element.atomic_number() > 0, "dummy element");
                numbers.push(element.atomic_number());
                let mut row = [0.0; 3];
                for coordinate in &mut row {
                    *coordinate = fields.next().context("missing coordinate")?.parse()?;
                    ensure!(f64::is_finite(*coordinate), "non-finite coordinate");
                }
                angstrom.push(row);
            }
            ensure!(
                numbers.len() == atoms && lines.all(|line| line.trim().is_empty()),
                "XYZ atom count does not match its coordinate block"
            );
            Ok((Self { numbers, angstrom }, comment))
        };
        read().with_context(|| format!("reading {}", path.display()))
    }

    pub fn write(&self, path: &Path, comment: &str) -> Result<()> {
        let mut xyz = format!("{}\n{comment} SPDX-License-Identifier: LGPL-3.0-or-later\n", self.numbers.len());
        for (&number, &[x, y, z]) in self.numbers.iter().zip(&self.angstrom) {
            let symbol = Element::from_atomic_number(number).context("invalid atomic number")?.symbol();
            writeln!(xyz, "{symbol} {x} {y} {z}")?;
        }
        fs::create_dir_all(parent(path))?;
        fs::write(path, xyz).with_context(|| format!("writing {}", path.display()))
    }

    pub fn bohr(&self) -> Vec<f64> {
        self.angstrom.as_flattened().iter().map(|x| x / ANGSTROM_PER_BOHR).collect()
    }

    pub fn with_bohr(&self, positions: &[f64]) -> Self {
        Self {
            numbers: self.numbers.clone(),
            angstrom: positions.as_chunks::<3>().0.iter().map(|row| row.map(|x| x * ANGSTROM_PER_BOHR)).collect(),
        }
    }
}

#[derive(Clone, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
pub struct Identity {
    pub case: String,
    pub set: String,
    pub start: String,
}

impl Identity {
    pub fn label(&self) -> String {
        format!("{}/{}.{}", self.set, self.case, self.start)
    }

    pub fn path(&self) -> PathBuf {
        PathBuf::from(self.label() + ".xyz")
    }
}

pub struct Start {
    pub id: Identity,
    pub geometry: Geometry,
    pub charge: i32,
    pub uhf: u32,
    pub smiles: Option<String>,
    pub atom_order: Option<String>,
    pub mmff_state: Option<String>,
}

impl Start {
    pub fn comment(&self) -> String {
        let mut text =
            format!("charge={} uhf={} smiles={}", self.charge, self.uhf, self.smiles.as_deref().unwrap_or(""));
        if let Some(order) = &self.atom_order {
            text.push_str(&format!(" atom_order={order}"));
        }
        if let Some(state) = &self.mmff_state {
            text.push_str(&format!(" mmff={state}"));
        }
        text
    }

    pub fn mmff_field(&self, geometry: &Geometry) -> Result<MolecularForceField> {
        ensure!(
            self.mmff_state.as_deref() != Some("unavailable"),
            "COSMolKit MMFF94 could not type this case during prepare"
        );
        let smiles = self.smiles.as_deref().context("MMFF94 needs a SMILES case")?;
        let molecule = molecule(smiles, self.charge)?;
        ensure!(
            numbers(&molecule) == geometry.numbers && self.atom_order.as_deref() == Some(&atom_order(&molecule)),
            "SMILES atom order differs from the order recorded by prepare"
        );
        let with_coordinates =
            molecule.with_only_3d_conformer(geometry.angstrom.iter().map(|row| row.to_vec()).collect())?;
        field(&with_coordinates)
    }
}

#[derive(clap::Args, Default)]
pub struct Selection {
    #[arg(long, value_delimiter = ',')]
    sets: Vec<String>,
    /// Case names or set/name pairs.
    #[arg(long, value_delimiter = ',')]
    cases: Vec<String>,
    #[arg(long, value_delimiter = ',', value_parser = ["etkdg", "mmff", "given", "rough"])]
    starts_kind: Vec<String>,
}

pub fn starts(root: &Path, selection: &Selection) -> Result<Vec<Start>> {
    let mut found = Vec::new();
    for set in fs::read_dir(root).with_context(|| format!("reading {}", root.display()))? {
        let set = set?;
        if !set.file_type()?.is_dir() {
            continue;
        }
        let name = set.file_name().into_string().map_err(|_| anyhow::anyhow!("non-UTF-8 set name"))?;
        if !selected(&selection.sets, &name) {
            continue;
        }
        validate_name(&name)?;
        for file in fs::read_dir(set.path())? {
            let file = file?;
            let path = file.path();
            if !file.file_type()?.is_file() || path.extension().and_then(|s| s.to_str()) != Some("xyz") {
                continue;
            }
            let (case, kind) = path
                .file_stem()
                .and_then(|s| s.to_str())
                .and_then(|s| s.rsplit_once('.'))
                .with_context(|| format!("expected <case>.<start>.xyz at {}", path.display()))?;
            if !(selection.cases.is_empty()
                || selection.cases.iter().any(|choice| choice == case || choice == &format!("{name}/{case}")))
                || !selected(&selection.starts_kind, kind)
            {
                continue;
            }
            ensure!(["etkdg", "mmff", "given", "rough"].contains(&kind), "unknown start kind {kind:?}");
            validate_name(case)?;
            let (geometry, comment) = Geometry::read(&path)?;
            let metadata = |key: &str| {
                comment
                    .split_whitespace()
                    .filter_map(|word| word.split_once('='))
                    .find_map(|(k, v)| (k == key && !v.is_empty()).then_some(v))
            };
            let charge = metadata("charge")
                .context("missing charge metadata")?
                .parse()
                .with_context(|| format!("invalid charge in {}", path.display()))?;
            let uhf = metadata("uhf")
                .context("missing uhf metadata")?
                .parse()
                .with_context(|| format!("invalid uhf in {}", path.display()))?;
            found.push(Start {
                id: Identity { case: case.into(), set: name.clone(), start: kind.into() },
                geometry,
                charge,
                uhf,
                smiles: metadata("smiles").map(str::to_owned),
                atom_order: metadata("atom_order").map(str::to_owned),
                mmff_state: metadata("mmff").map(str::to_owned),
            });
        }
    }
    found.sort_by(|a, b| a.id.cmp(&b.id));
    ensure!(!found.is_empty(), "no starts matched in {}", root.display());
    Ok(found)
}

fn selected(values: &[String], item: &str) -> bool {
    values.is_empty() || values.iter().any(|value| value == item)
}

fn validate_name(name: &str) -> Result<()> {
    ensure!(
        !name.is_empty() && name.bytes().all(|c| c.is_ascii_alphanumeric() || b"_-".contains(&c)),
        "invalid case or set name {name:?}"
    );
    Ok(())
}

pub fn parent(path: &Path) -> &Path {
    path.parent().filter(|parent| !parent.as_os_str().is_empty()).unwrap_or_else(|| Path::new("."))
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Manifest {
    cases: Vec<Case>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Case {
    name: String,
    set: String,
    smiles: Option<String>,
    xyz: Option<PathBuf>,
    charge: i32,
    uhf: u32,
}

fn molecule(smiles: &str, charge: i32) -> Result<Molecule> {
    let molecule = Molecule::from_smiles(smiles)?.with_hydrogens()?.sanitize()?;
    let total: i32 = molecule.atoms().iter().map(|atom| i32::from(atom.formal_charge())).sum();
    ensure!(charge == total, "manifest charge {charge} differs from SMILES charge {total}");
    Ok(molecule)
}

fn numbers(molecule: &Molecule) -> Vec<u8> {
    molecule.atoms().iter().map(|atom| atom.atomic_number()).collect()
}

// Record indexed connectivity in addition to the element sequence, so reordering unlike
// environments of the same element is detectable. FNV-1a and the text encoding are fixed.
fn atom_order(molecule: &Molecule) -> String {
    let mut topology = format!("{:?}|", numbers(molecule));
    let mut bonds = molecule
        .bonds()
        .iter()
        .map(|bond| {
            let a = bond.begin().index();
            let b = bond.end().index();
            (a.min(b), a.max(b), format!("{:?}", bond.order()))
        })
        .collect::<Vec<_>>();
    bonds.sort_unstable();
    for (a, b, order) in bonds {
        write!(topology, "{a}-{b}:{order};").expect("writing to String");
    }
    format!("{:016x}", stable_hash(topology.as_bytes()))
}

fn field(molecule: &Molecule) -> Result<MolecularForceField> {
    let params = MmffForceFieldParams::new(None, "MMFF94".into(), f64::INFINITY, true);
    Ok(molecule.mmff_force_field_with_params(&params)?)
}

pub fn prepare(manifest: &Path, out: &Path) -> Result<()> {
    let cases: Manifest =
        toml::from_str(&fs::read_to_string(manifest)?).with_context(|| format!("reading {}", manifest.display()))?;
    let mut names = HashSet::new();
    for case in &cases.cases {
        validate_name(&case.set)?;
        validate_name(&case.name)?;
        ensure!(names.insert((&case.set, &case.name)), "duplicate case {}/{}", case.set, case.name);
        ensure!(
            case.smiles.is_some() != case.xyz.is_some(),
            "{}/{} needs exactly one of smiles or xyz",
            case.set,
            case.name
        );
        if let Some(smiles) = &case.smiles {
            ensure!(!smiles.is_empty() && !smiles.contains(char::is_whitespace), "SMILES must be a single token");
        }
    }
    let mut counts = std::collections::BTreeMap::<&str, [usize; 4]>::new();
    for case in &cases.cases {
        let count = counts.entry(&case.set).or_default();
        count[0] += 1;
        let (mut first, second) = if let Some(smiles) = &case.smiles {
            let prepared = match embed(smiles, case.charge) {
                Ok(prepared) => prepared,
                Err(error) => {
                    eprintln!("Skipped {}/{}: {error:#}", case.set, case.name);
                    continue;
                }
            };
            let (geometry, atom_order, embedded) = prepared;
            let (relaxed, state) = match field(&embedded) {
                Ok(mut force_field) => match minimize(&mut force_field) {
                    Ok(true) => (force_field.positions(), "minimized"),
                    Ok(false) => {
                        eprintln!(
                            "{}/{}: MMFF94 minimizer reached its limit; retaining its last geometry",
                            case.set, case.name
                        );
                        count[3] += 1;
                        (force_field.positions(), "max_iterations")
                    }
                    Err(error) => {
                        eprintln!(
                            "{}/{}: MMFF94 minimization failed; using ETKDG for mmff start: {error:#}",
                            case.set, case.name
                        );
                        count[3] += 1;
                        (geometry.angstrom.clone(), "minimization_failed")
                    }
                },
                Err(error) => {
                    eprintln!("{}/{}: MMFF94 unavailable; using ETKDG for mmff start: {error:#}", case.set, case.name);
                    count[3] += 1;
                    (geometry.angstrom.clone(), "unavailable")
                }
            };
            let start = Start {
                id: Identity { case: case.name.clone(), set: case.set.clone(), start: "etkdg".into() },
                geometry,
                charge: case.charge,
                uhf: case.uhf,
                smiles: Some(smiles.clone()),
                atom_order: Some(atom_order),
                mmff_state: Some(state.into()),
            };
            (start, ("mmff", relaxed))
        } else {
            let Some(path) = &case.xyz else {
                bail!("missing XYZ path");
            };
            let (geometry, _) = Geometry::read(&parent(manifest).join(path))?;
            let rough = perturb(&geometry.angstrom, &case.name);
            let start = Start {
                id: Identity { case: case.name.clone(), set: case.set.clone(), start: "given".into() },
                geometry,
                charge: case.charge,
                uhf: case.uhf,
                smiles: None,
                atom_order: None,
                mmff_state: None,
            };
            (start, ("rough", rough))
        };
        first.geometry.write(&out.join(first.id.path()), &first.comment())?;
        first.id.start = second.0.into();
        first.geometry.angstrom = second.1;
        first.geometry.write(&out.join(first.id.path()), &first.comment())?;
        count[1] += 1;
        count[2] += 2;
    }
    println!("set\tcases\tprepared\tstarts\tmmff_warnings");
    for (set, [cases, prepared, starts, mmff_warnings]) in counts {
        println!("{set}\t{cases}\t{prepared}\t{starts}\t{mmff_warnings}");
    }
    Ok(())
}

fn embed(smiles: &str, charge: i32) -> Result<(Geometry, String, Molecule)> {
    let molecule = molecule(smiles, charge)?;
    let order = atom_order(&molecule);
    let mut params = EmbedParams::etkdg_v3();
    params.random_seed = 42;
    params.num_threads = 1;
    let embedded = molecule.with_3d_conformer_result_with_params(&params)?.molecule;
    ensure!(order == atom_order(&embedded), "ETKDG changed atom order");
    let [conformer] = embedded.conformers_3d() else {
        bail!("ETKDG did not return exactly one conformer");
    };
    let geometry = Geometry { numbers: numbers(&embedded), angstrom: conformer.coordinates().to_vec() };
    Ok((geometry, order, embedded))
}

fn minimize(field: &mut MolecularForceField) -> Result<bool> {
    let defaults = ForceFieldMinimizeParams::default();
    let params = ForceFieldMinimizeParams::new(10_000, defaults.force_tolerance(), defaults.energy_tolerance());
    Ok(field.minimize_with_params_(&params)?.converged())
}

fn stable_hash(input: &[u8]) -> u64 {
    input.iter().fold(0xcbf29ce484222325, |hash, &byte| (hash ^ u64::from(byte)).wrapping_mul(0x100000001b3))
}

fn perturb(coordinates: &[[f64; 3]], name: &str) -> Vec<[f64; 3]> {
    let mut state = stable_hash(name.as_bytes());
    let mut uniform = || {
        // SplitMix64, with a uniform variate in (0, 1] for the Box-Muller transform.
        state = state.wrapping_add(0x9e3779b97f4a7c15);
        let mut bits = state;
        bits = (bits ^ (bits >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
        bits = (bits ^ (bits >> 27)).wrapping_mul(0x94d049bb133111eb);
        bits ^= bits >> 31;
        ((bits >> 11) + 1) as f64 / (1u64 << 53) as f64
    };
    coordinates
        .iter()
        .map(|row| {
            row.map(|x| {
                let radius = (-2.0 * uniform().ln()).sqrt();
                x + 0.10 * radius * (std::f64::consts::TAU * uniform()).cos()
            })
        })
        .collect()
}
