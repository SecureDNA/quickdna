#![feature(array_chunks)]
#![allow(non_snake_case)]

use pyo3::{exceptions::PyValueError, prelude::*, types::PyBytes};
use trans_table::TranslationTable;

pub mod trans_table; // needs to be public for bin/gen_table

impl std::convert::From<trans_table::TranslationError> for PyErr {
    fn from(err: trans_table::TranslationError) -> PyErr {
        PyValueError::new_err(err.to_string())
    }
}

#[pyfunction]
fn _check_table(table: u8) -> PyResult<()> {
    let _ = TranslationTable::try_from(table)?;
    Ok(())
}

#[pyfunction]
fn _translate(py: Python, table: u8, dna: &PyBytes) -> PyResult<PyObject> {
    let bytes = trans_table::translate(table, dna.as_bytes())?;
    Ok(PyBytes::new(py, &bytes).into())
}

#[pyfunction]
fn _reverse_complement(py: Python, dna: &PyBytes) -> PyResult<PyObject> {
    let bytes = trans_table::reverse_complement(dna.as_bytes())?;
    Ok(PyBytes::new(py, &bytes).into())
}

#[pymodule]
fn quickdna(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(_check_table, m)?)?;
    m.add_function(wrap_pyfunction!(_translate, m)?)?;
    m.add_function(wrap_pyfunction!(_reverse_complement, m)?)?;

    Ok(())
}
