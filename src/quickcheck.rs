// Copyright 2021-2024 SecureDNA Stiftung (SecureDNA Foundation) <licensing@securedna.org>
// SPDX-License-Identifier: MIT OR Apache-2.0

use quickcheck::{Arbitrary, Gen};

use crate::{Codon, CodonAmbiguous, DnaSequence, Nucleotide, NucleotideAmbiguous, NucleotideLike};

impl Arbitrary for Nucleotide {
    fn arbitrary(g: &mut Gen) -> Self {
        *g.choose(&[Self::A, Self::T, Self::C, Self::G])
            .expect("Gen should be able to choose a Nucleotide")
    }
}

impl Arbitrary for NucleotideAmbiguous {
    fn arbitrary(g: &mut Gen) -> Self {
        *g.choose(&[
            Self::A,
            Self::T,
            Self::C,
            Self::G,
            Self::W,
            Self::M,
            Self::R,
            Self::Y,
            Self::S,
            Self::K,
            Self::B,
            Self::V,
            Self::D,
            Self::H,
            Self::N,
        ])
        .expect("Gen should be able to choose a NucleotideAmbiguous")
    }
}

impl Arbitrary for Codon {
    fn arbitrary(g: &mut Gen) -> Self {
        Self([
            Nucleotide::arbitrary(g),
            Nucleotide::arbitrary(g),
            Nucleotide::arbitrary(g),
        ])
    }
}

impl Arbitrary for CodonAmbiguous {
    fn arbitrary(g: &mut Gen) -> Self {
        Self([
            NucleotideAmbiguous::arbitrary(g),
            NucleotideAmbiguous::arbitrary(g),
            NucleotideAmbiguous::arbitrary(g),
        ])
    }
}

impl<T: Arbitrary + NucleotideLike> Arbitrary for DnaSequence<T> {
    fn arbitrary(g: &mut Gen) -> Self {
        Self::new(Arbitrary::arbitrary(g))
    }
}
