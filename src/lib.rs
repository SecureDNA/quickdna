#![feature(array_chunks)]

use pyo3::{exceptions::PyValueError, prelude::*, types::PyBytes};
use trans_table::TranslationTable;

pub mod trans_table; // needs to be public for bin/gen_table

impl std::convert::From<trans_table::TranslationError> for PyErr {
    fn from(err: trans_table::TranslationError) -> PyErr {
        PyValueError::new_err(err.to_string())
    }
}

#[pyfunction]
fn check_table(table: u8) -> PyResult<()> {
    let _ = TranslationTable::try_from(table)?;
    Ok(())
}

#[pyfunction]
fn translate(py: Python, table: u8, dna: &PyBytes) -> PyResult<PyObject> {
    let bytes = trans_table::translate(table, dna.as_bytes())?;
    Ok(PyBytes::new(py, &bytes).into())
}

#[pyfunction]
fn reverse_complement(py: Python, dna: &PyBytes) -> PyResult<PyObject> {
    let bytes = trans_table::reverse_complement(dna.as_bytes())?;
    Ok(PyBytes::new(py, &bytes).into())
}

#[pymodule]
fn _quickdna(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(check_table, m)?)?;
    m.add_function(wrap_pyfunction!(translate, m)?)?;
    m.add_function(wrap_pyfunction!(reverse_complement, m)?)?;

    Ok(())
}
