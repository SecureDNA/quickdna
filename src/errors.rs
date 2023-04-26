// Copyright 2021-2023 SecureDNA Stiftung (SecureDNA Foundation) <license@securedna.org>
// SPDX-License-Identifier: MIT OR Apache-2.0

use std::array::TryFromSliceError;

use thiserror::Error;

#[derive(Debug, Clone, Error)]
#[error("on line {line_number}: {error}")]
pub struct Located<E> {
    pub line_number: usize,

    #[source]
    pub error: E,
}

#[derive(Debug, Clone, Error)]
pub enum TranslationError {
    #[error("non-ascii byte: {:x?}", .0)]
    NonAsciiByte(u8),
    #[error("non-ascii char: {:x?}", .0)]
    NonAsciiChar(char),
    #[error("bad nucleotide: {:?}", .0)]
    BadNucleotide(char),
    #[error("unexpected ambiguous nucleotide: {:?}", .0)]
    UnexpectedAmbiguousNucleotide(char),
    #[error("not a ncbi translation table: {}", .0)]
    BadTranslationTable(u8),
}

#[derive(Debug, Clone, Error)]
pub enum CodonError {
    #[error("{:?}", .0)]
    BadTranslation(#[from] TranslationError),
    #[error("{:?}", .0)]
    BadSlice(#[from] TryFromSliceError),
}
