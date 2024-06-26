// Copyright 2021-2024 SecureDNA Stiftung (SecureDNA Foundation) <licensing@securedna.org>
// SPDX-License-Identifier: MIT OR Apache-2.0

//! This binary is used to generate the src/tables.dat datafile.
//! Should be run from the repository root.
//! Translation table via Wikipedia.

use std::{collections::HashSet, fs};

use quickdna::{
    trans_table::{CodonIdx, TranslationTable},
    CodonAmbiguous, Nucleotide, NucleotideAmbiguous, NucleotideLike,
};

// via https://en.wikipedia.org/wiki/List_of_genetic_codes
const TABLE: &str = "
ATT || I || I || I || I || I || I || I || I || I || I || I || I || I || I || I || I || I || I || I || I || I || I || I || I || I || I || I
ATC || I || I || I || I || I || I || I || I || I || I || I || I || I || I || I || I || I || I || I || I || I || I || I || I || I || I || I
ATA || I || M || M || I || M || I || I || I || I || I || M || I || I || I || M || I || I || I || I || I || I || I || I || I || I || I || I
ATG || M || M || M || M || M || M || M || M || M || M || M || M || M || M || M || M || M || M || M || M || M || M || M || M || M || M || M
ACT || T || T || T || T || T || T || T || T || T || T || T || T || T || T || T || T || T || T || T || T || T || T || T || T || T || T || T
ACC || T || T || T || T || T || T || T || T || T || T || T || T || T || T || T || T || T || T || T || T || T || T || T || T || T || T || T
ACA || T || T || T || T || T || T || T || T || T || T || T || T || T || T || T || T || T || T || T || T || T || T || T || T || T || T || T
ACG || T || T || T || T || T || T || T || T || T || T || T || T || T || T || T || T || T || T || T || T || T || T || T || T || T || T || T
AAT || N || N || N || N || N || N || N || N || N || N || N || N || N || N || N || N || N || N || N || N || N || N || N || N || N || N || N
AAC || N || N || N || N || N || N || N || N || N || N || N || N || N || N || N || N || N || N || N || N || N || N || N || N || N || N || N
AAA || K || K || K || K || K || K || N || K || K || K || K || N || K || K || N || K || K || K || K || K || K || K || K || K || K || K || K
AAG || K || K || K || K || K || K || K || K || K || K || K || K || K || K || K || K || K || K || K || K || K || K || K || K || K || K || K
AGT || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S
AGC || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S
AGA || R || * || R || R || S || R || S || R || R || R || G || S || R || R || S || R || R || S || R || R || R || R || R || R || R || R || S
AGG || R || * || R || R || S || R || S || R || R || R || G || S || R || R || S || R || R || K || R || R || R || R || R || R || R || R || K
TAT || Y || Y || Y || Y || Y || Y || Y || Y || Y || Y || Y || Y || Y || Y || Y || Y || Y || Y || Y || Y || Y || Y || Y || Y || Y || Y || Y
TAC || Y || Y || Y || Y || Y || Y || Y || Y || Y || Y || Y || Y || Y || Y || Y || Y || Y || Y || Y || Y || Y || Y || Y || Y || Y || Y || Y
TAA || * || * || * || * || * || Q || * || * || * || * || * || Y || * || * || * || * || * || * || * || * || Q || Q || Y || E || E || * || Y
TAG || * || * || * || * || * || Q || * || * || * || * || * || * || Q || L || * || L || * || * || * || * || Q || Q || Y || E || E || W || *
TTA || L || L || L || L || L || L || L || L || L || L || L || L || L || L || L || L || * || L || L || L || L || L || L || L || L || L || L
TTT || F || F || F || F || F || F || F || F || F || F || F || F || F || F || F || F || F || F || F || F || F || F || F || F || F || F || F
TTC || F || F || F || F || F || F || F || F || F || F || F || F || F || F || F || F || F || F || F || F || F || F || F || F || F || F || F
TTG || L || L || L || L || L || L || L || L || L || L || L || L || L || L || L || L || L || L || L || L || L || L || L || L || L || L || L
TCT || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S
TCC || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S
TCA || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || * || S || S || S || S || S || S || S || S || S || S || S
TCG || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S || S
TGT || C || C || C || C || C || C || C || C || C || C || C || C || C || C || C || C || C || C || C || C || C || C || C || C || C || C || C
TGC || C || C || C || C || C || C || C || C || C || C || C || C || C || C || C || C || C || C || C || C || C || C || C || C || C || C || C
TGA || * || W || W || W || W || * || W || C || * || * || W || W || * || * || W || * || * || W || G || * || W || W || * || * || W || * || W
TGG || W || W || W || W || W || W || W || W || W || W || W || W || W || W || W || W || W || W || W || W || W || W || W || W || W || W || W
CTT || L || L || T || L || L || L || L || L || L || L || L || L || L || L || L || L || L || L || L || L || L || L || L || L || L || L || L
CTC || L || L || T || L || L || L || L || L || L || L || L || L || L || L || L || L || L || L || L || L || L || L || L || L || L || L || L
CTA || L || L || T || L || L || L || L || L || L || L || L || L || L || L || L || L || L || L || L || L || L || L || L || L || L || L || L
CTG || L || L || T || L || L || L || L || L || L || S || L || L || L || L || L || L || L || L || L || A || L || L || L || L || L || L || L
CCT || P || P || P || P || P || P || P || P || P || P || P || P || P || P || P || P || P || P || P || P || P || P || P || P || P || P || P
CCC || P || P || P || P || P || P || P || P || P || P || P || P || P || P || P || P || P || P || P || P || P || P || P || P || P || P || P
CCA || P || P || P || P || P || P || P || P || P || P || P || P || P || P || P || P || P || P || P || P || P || P || P || P || P || P || P
CCG || P || P || P || P || P || P || P || P || P || P || P || P || P || P || P || P || P || P || P || P || P || P || P || P || P || P || P
CAT || H || H || H || H || H || H || H || H || H || H || H || H || H || H || H || H || H || H || H || H || H || H || H || H || H || H || H
CAC || H || H || H || H || H || H || H || H || H || H || H || H || H || H || H || H || H || H || H || H || H || H || H || H || H || H || H
CAA || Q || Q || Q || Q || Q || Q || Q || Q || Q || Q || Q || Q || Q || Q || Q || Q || Q || Q || Q || Q || Q || Q || Q || Q || Q || Q || Q
CAG || Q || Q || Q || Q || Q || Q || Q || Q || Q || Q || Q || Q || Q || Q || Q || Q || Q || Q || Q || Q || Q || Q || Q || Q || Q || Q || Q
CGT || R || R || R || R || R || R || R || R || R || R || R || R || R || R || R || R || R || R || R || R || R || R || R || R || R || R || R
CGC || R || R || R || R || R || R || R || R || R || R || R || R || R || R || R || R || R || R || R || R || R || R || R || R || R || R || R
CGA || R || R || R || R || R || R || R || R || R || R || R || R || R || R || R || R || R || R || R || R || R || R || R || R || R || R || R
CGG || R || R || R || R || R || R || R || R || R || R || R || R || R || R || R || R || R || R || R || R || R || R || R || R || R || R || R
GTT || V || V || V || V || V || V || V || V || V || V || V || V || V || V || V || V || V || V || V || V || V || V || V || V || V || V || V
GTC || V || V || V || V || V || V || V || V || V || V || V || V || V || V || V || V || V || V || V || V || V || V || V || V || V || V || V
GTA || V || V || V || V || V || V || V || V || V || V || V || V || V || V || V || V || V || V || V || V || V || V || V || V || V || V || V
GTG || V || V || V || V || V || V || V || V || V || V || V || V || V || V || V || V || V || V || V || V || V || V || V || V || V || V || V
GCT || A || A || A || A || A || A || A || A || A || A || A || A || A || A || A || A || A || A || A || A || A || A || A || A || A || A || A
GCC || A || A || A || A || A || A || A || A || A || A || A || A || A || A || A || A || A || A || A || A || A || A || A || A || A || A || A
GCA || A || A || A || A || A || A || A || A || A || A || A || A || A || A || A || A || A || A || A || A || A || A || A || A || A || A || A
GCG || A || A || A || A || A || A || A || A || A || A || A || A || A || A || A || A || A || A || A || A || A || A || A || A || A || A || A
GAT || D || D || D || D || D || D || D || D || D || D || D || D || D || D || D || D || D || D || D || D || D || D || D || D || D || D || D
GAC || D || D || D || D || D || D || D || D || D || D || D || D || D || D || D || D || D || D || D || D || D || D || D || D || D || D || D
GAA || E || E || E || E || E || E || E || E || E || E || E || E || E || E || E || E || E || E || E || E || E || E || E || E || E || E || E
GAG || E || E || E || E || E || E || E || E || E || E || E || E || E || E || E || E || E || E || E || E || E || E || E || E || E || E || E
GGT || G || G || G || G || G || G || G || G || G || G || G || G || G || G || G || G || G || G || G || G || G || G || G || G || G || G || G
GGC || G || G || G || G || G || G || G || G || G || G || G || G || G || G || G || G || G || G || G || G || G || G || G || G || G || G || G
GGA || G || G || G || G || G || G || G || G || G || G || G || G || G || G || G || G || G || G || G || G || G || G || G || G || G || G || G
GGG || G || G || G || G || G || G || G || G || G || G || G || G || G || G || G || G || G || G || G || G || G || G || G || G || G || G || G
";

fn ambiguous_codon_protein(
    codon: CodonAmbiguous,
    table_idx: usize,
    translation_tables: &[u8; TranslationTable::LOOKUP_SIZE],
) -> u8 {
    let mut seen_proteins: HashSet<u8> = HashSet::new();

    for &n1 in codon.0[0].possibilities() {
        for &n2 in codon.0[1].possibilities() {
            for &n3 in codon.0[2].possibilities() {
                let codon_idx: usize = CodonIdx::from([n1, n2, n3]).into();
                let candidate_protein =
                    translation_tables[table_idx * TranslationTable::CODONS_PER_TABLE + codon_idx];

                seen_proteins.insert(candidate_protein);
                if seen_proteins.len() >= 3 {
                    return b'X';
                }
            }
        }
    }

    let mut seen_vec: Vec<u8> = seen_proteins.into_iter().collect();
    seen_vec.sort();

    match seen_vec.as_slice() {
        [single] => *single,
        [b'D', b'N'] => b'B', // B = Asx = Asparagine or Aspartic acid
        [b'E', b'Q'] => b'Z', // Z = Glx = Glutamine or Glutamic acid
        [b'I', b'L'] => b'J', // J = Xle = Leucine or Isoleucine
        _ => b'X',
    }
}

fn gen_translation_tables() -> [u8; TranslationTable::LOOKUP_SIZE] {
    let mut translation_tables = [b'*'; TranslationTable::LOOKUP_SIZE];

    for line in TABLE.trim().split('\n') {
        let mut cols = line.split(" || ");

        let codon = cols.next().unwrap().trim().as_bytes();
        let codon_idx: usize = CodonIdx::from([
            Nucleotide::try_from(codon[0]).expect("failed to parse nucleotide 1 in codon"),
            Nucleotide::try_from(codon[1]).expect("failed to parse nucleotide 2 in codon"),
            Nucleotide::try_from(codon[2]).expect("failed to parse nucleotide 3 in codon"),
        ])
        .into();

        for (table_idx, aa) in cols.enumerate() {
            let aa = aa.trim().as_bytes()[0];
            translation_tables[table_idx * TranslationTable::CODONS_PER_TABLE + codon_idx] = aa;
        }
    }

    // setup ambiguous codons
    for a in NucleotideAmbiguous::ALL {
        for b in NucleotideAmbiguous::ALL {
            for c in NucleotideAmbiguous::ALL {
                if a.is_ambiguous() || b.is_ambiguous() || c.is_ambiguous() {
                    let codon = CodonAmbiguous([a, b, c]);
                    let codon_idx: usize = CodonIdx::from(codon).into();

                    for table_idx in 0..TranslationTable::N_TRANS_TABLES {
                        let protein =
                            ambiguous_codon_protein(codon, table_idx, &translation_tables);
                        translation_tables
                            [table_idx * TranslationTable::CODONS_PER_TABLE + codon_idx] = protein;
                    }
                }
            }
        }
    }

    translation_tables
}

fn main() {
    let tables = gen_translation_tables();
    fs::write("src/tables.dat", tables).expect("failed to write tables.dat");
    println!("Wrote codon translation table data to src/tables.dat.");
}
