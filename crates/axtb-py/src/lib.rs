// SPDX-License-Identifier: LGPL-3.0-or-later
//! The extension module `axtb._axtb`: the private engine under axtb's Python interface
//! (`axtb.compute`, `axtb.Session`, `axtb.program`; python/axtb).
//!
//! The Python layer checks a `qcdata.ProgramInput` and turns it into the driver's flat request;
//! this module runs it on `axtb-driver` (`engine`: one calculation as the CLI runs it, and a
//! session's calculator with one committed wavefunction per lineage) and returns the response as
//! a dict of numpy arrays. It also answers the planning's questions from the driver's tables
//! (canonical names, the elements and shells of a method, numerics identities, the one-thread
//! sizes). Errors are `axtb._axtb.Error` with the driver's kind and code.

mod engine;

use pyo3::prelude::*;
use pyo3::types::{PyBool, PyDict, PyList};

/// The JSON value as Python objects (null, bool, number, string, list, dict).
pub(crate) fn to_py(py: Python<'_>, v: &serde_json::Value) -> PyResult<PyObject> {
    use serde_json::Value;
    Ok(match v {
        Value::Null => py.None(),
        Value::Bool(b) => PyBool::new(py, *b).to_owned().into_any().unbind(),
        Value::Number(n) => match (n.as_u64(), n.as_i64(), n.as_f64()) {
            (Some(u), _, _) => u.into_pyobject(py)?.into_any().unbind(),
            (_, Some(i), _) => i.into_pyobject(py)?.into_any().unbind(),
            (_, _, Some(f)) => f.into_pyobject(py)?.into_any().unbind(),
            _ => py.None(),
        },
        Value::String(s) => s.into_pyobject(py)?.into_any().unbind(),
        Value::Array(a) => {
            let l = PyList::empty(py);
            for x in a {
                l.append(to_py(py, x)?)?;
            }
            l.into_any().unbind()
        }
        Value::Object(m) => {
            let d = PyDict::new(py);
            for (k, x) in m {
                d.set_item(k, to_py(py, x)?)?;
            }
            d.into_any().unbind()
        }
    })
}

/// axtb's extension module. `axtb/__init__.py` checks the CPU before it imports this module (an
/// x86-64-v3 build runs AVX instructions as soon as the module loads).
#[pymodule]
fn _axtb(m: &Bound<'_, PyModule>) -> PyResult<()> {
    let py = m.py();
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add("Error", py.get_type::<engine::Error>())?;
    m.add("ANALYSIS_NUMERICS", axtb_driver::ANALYSIS_NUMERICS)?;
    m.add_class::<engine::Engine>()?;
    m.add_function(wrap_pyfunction!(engine::compute_once, m)?)?;
    m.add_function(wrap_pyfunction!(engine::canonical, m)?)?;
    m.add_function(wrap_pyfunction!(engine::supported_elements, m)?)?;
    m.add_function(wrap_pyfunction!(engine::numerics, m)?)?;
    m.add_function(wrap_pyfunction!(engine::single_thread_atoms, m)?)?;
    m.add_function(wrap_pyfunction!(engine::valence_electrons, m)?)?;
    m.add_function(wrap_pyfunction!(engine::shell_layout, m)?)?;
    m.add_function(wrap_pyfunction!(engine::solvents, m)?)?;
    m.add_function(wrap_pyfunction!(engine::canonical_solvent, m)?)?;
    m.add_function(wrap_pyfunction!(engine::lapack_warning, m)?)?;
    Ok(())
}
