// Copyright 2021-2024 SecureDNA Stiftung (SecureDNA Foundation) <licensing@securedna.org>
// SPDX-License-Identifier: MIT OR Apache-2.0

/// A trait used by the FASTA parser: `T: Extendable` is the type of the
/// contents of the FASTA file (String or DnaSequence or ProteinSequence).
/// This trait lets us generically concatenate parsed content lines.
pub trait Extendable: Default {
    /// Is this value empty when ignoring whitespace?
    fn is_blank(&self) -> bool;

    /// Extend this value by another value.
    fn extend(&mut self, other: Self);
}

impl Extendable for String {
    fn is_blank(&self) -> bool {
        self.trim().is_empty()
    }

    fn extend(&mut self, other: Self) {
        self.push_str(&other)
    }
}
