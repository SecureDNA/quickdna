// Copyright 2021-2023 SecureDNA Stiftung (SecureDNA Foundation) <licensing@securedna.org>
// SPDX-License-Identifier: MIT OR Apache-2.0

use quickdna::{Nucleotide, NucleotideAmbiguous, NucleotideIter, NucleotideLike, TranslationTable};

// The purpose of this is to ensure that it's possible to write code that's
// generic on NucleotideLike types, without excessive where-clause constraints.
// In particular:
// * This only relies on being given a slice of NucleotideLikes
// * We're able to use NucleotideIter (though it sadly requires .copied())
// * We're able to map the resulting generic Codons through a translation table.
fn get_rc_peptide_from_generic_slice(dna: &[impl NucleotideLike]) -> String {
    let table = TranslationTable::Ncbi1.to_fn();
    let peptide = dna
        .iter()
        .copied()
        .reverse_complement()
        .codons()
        .map(table)
        .collect();
    String::from_utf8(peptide).unwrap()
}

#[test]
fn test_generic_method() {
    let dna = {
        use Nucleotide::*;
        [C, A, T, T, A, G]
    };
    let aas = get_rc_peptide_from_generic_slice(&dna);
    assert_eq!(aas, "LM");

    let dna_ambiguous = {
        use NucleotideAmbiguous::*;
        [C, A, T, T, A, G]
    };
    let aas_ambiguous = get_rc_peptide_from_generic_slice(&dna_ambiguous);
    assert_eq!(aas_ambiguous, "LM");
}
