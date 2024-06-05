#![allow(clippy::borrow_deref_ref)] // TODO: broken clippy lint?
                                    // Copyright 2021-2024 SecureDNA Stiftung (SecureDNA Foundation) <licensing@securedna.org>
                                    // SPDX-License-Identifier: MIT OR Apache-2.0

use pyo3::{exceptions::PyValueError, prelude::*, types::PyBytes};

use crate::{
    errors::TranslationError,
    trans_table::{reverse_complement_bytes, TranslationTable},
    Nucleotide, NucleotideAmbiguous,
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

/// Translate a bytestring of DNA nucleotides into a bytestring of amino acids.
///
/// The input string is allowed to contain IUPAC ambiguity codes; ambiguous
/// codons are represented by `X` in the output.
///
/// * `translate(b"CCNTACACK CATNCNAAT")` returns `b"PYTHXN"`
#[pyfunction]
fn _translate(py: Python, table: u8, dna: &PyBytes) -> PyResult<PyObject> {
    let table = TranslationTable::try_from(table)?;
    let bytes = table.translate_dna_bytes::<NucleotideAmbiguous>(dna.as_bytes())?;
    Ok(PyBytes::new(py, &bytes).into())
}

/// Translate a bytestring of DNA nucleotides into a bytestring of amino acids.
///
/// The input string is validated to consist of unambiguous nucleotides (no IUPAC ambiguity codes).
///
/// * `translate_strict(b"AAACCCTTTGGG")` returns `b"KPFG"`
/// * `translate_strict(b"AAACCCTTTGGN")` is an error.
#[pyfunction]
fn _translate_strict(py: Python, table: u8, dna: &PyBytes) -> PyResult<PyObject> {
    let table = TranslationTable::try_from(table)?;
    let bytes = table.translate_dna_bytes::<Nucleotide>(dna.as_bytes())?;
    Ok(PyBytes::new(py, &bytes).into())
}

/// Get the reverse complement of a bytestring of DNA nucleotides.
///
/// The input string is allowed to contain IUPAC ambiguity codes.
///
/// * `reverse_complement(b"AAAAABCCC")` returns `b"GGGVTTTTT"`
#[pyfunction]
fn _reverse_complement(py: Python, dna: &PyBytes) -> PyResult<PyObject> {
    let bytes = reverse_complement_bytes::<NucleotideAmbiguous>(dna.as_bytes())?;
    Ok(PyBytes::new(py, &bytes).into())
}

/// Get the reverse complement of a bytestring of DNA nucleotides.
///
/// The input string is validated to consist of unambiguous nucleotides (no IUPAC ambiguity codes).
///
/// * `reverse_complement_strict(b"AAAAAACCC")` returns `b"GGGTTTTTT"`
/// * `reverse_complement_strict(b"AAAAAACCN")` is an error.
#[pyfunction]
fn _reverse_complement_strict(py: Python, dna: &PyBytes) -> PyResult<PyObject> {
    let bytes = reverse_complement_bytes::<Nucleotide>(dna.as_bytes())?;
    Ok(PyBytes::new(py, &bytes).into())
}

#[pymodule]
fn quickdna(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(_check_table, m)?)?;
    m.add_function(wrap_pyfunction!(_translate, m)?)?;
    m.add_function(wrap_pyfunction!(_translate_strict, m)?)?;
    m.add_function(wrap_pyfunction!(_reverse_complement, m)?)?;
    m.add_function(wrap_pyfunction!(_reverse_complement_strict, m)?)?;

    Ok(())
}
