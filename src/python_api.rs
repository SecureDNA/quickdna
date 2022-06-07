#![allow(clippy::borrow_deref_ref)]

use pyo3::{exceptions::PyValueError, prelude::*, types::PyBytes};

use crate::{
    errors::TranslationError,
    trans_table::{reverse_complement_bytes, TranslationTable},
};

impl From<TranslationError> for PyErr {
    fn from(err: TranslationError) -> PyErr {
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
    let table = TranslationTable::try_from(table)?;
    let bytes = table.translate_dna_bytes(dna.as_bytes())?;
    Ok(PyBytes::new(py, &bytes).into())
}

#[pyfunction]
fn _reverse_complement(py: Python, dna: &PyBytes) -> PyResult<PyObject> {
    let bytes = reverse_complement_bytes(dna.as_bytes())?;
    Ok(PyBytes::new(py, &bytes).into())
}

#[pymodule]
fn quickdna(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(_check_table, m)?)?;
    m.add_function(wrap_pyfunction!(_translate, m)?)?;
    m.add_function(wrap_pyfunction!(_reverse_complement, m)?)?;

    Ok(())
}
