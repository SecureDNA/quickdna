//! This binary is used to generate the src/tables.dat datafile.
//! Should be run from the repository root.
//! Translation table via Wikipedia.

use std::fs;

use quickdna::{
    trans_table::{CodonIdx, TranslationTable},
    Codon, Nucleotide,
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

fn codon_first_n_index(codon: Codon) -> Option<usize> {
    let codon: [Nucleotide; 3] = codon.into();
    codon
        .iter()
        .enumerate()
        .find_map(|(idx, &n)| if n == Nucleotide::N { Some(idx) } else { None })
}

fn ambiguous_codon_protein(
    codon: Codon,
    table_idx: usize,
    translation_tables: [u8; TranslationTable::LOOKUP_SIZE],
) -> u8 {
    match codon.count(Nucleotide::N) {
        0 => panic!("not an ambigious codon: {}", codon),
        1 => {
            // try and see if all possible permutations of this codon map to the same protein
            let mut seen_protein = None;
            let n_index = codon_first_n_index(codon).unwrap();

            for possible in [Nucleotide::A, Nucleotide::T, Nucleotide::C, Nucleotide::G] {
                let mut candidate: [Nucleotide; 3] = codon.into();
                candidate[n_index] = possible;
                let codon_idx: usize = CodonIdx::from(candidate).into();
                let candidate_protein =
                    translation_tables[table_idx * TranslationTable::CODONS_PER_TABLE + codon_idx];

                match seen_protein {
                    None => {
                        seen_protein = Some(candidate_protein);
                    }
                    Some(p) if p == candidate_protein => {}
                    _ => return b'X',
                }
            }

            seen_protein.unwrap()
        }
        // all codons with 2 or more Ns always map to X
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
    for &a in Nucleotide::NUCLEOTIDES.iter() {
        for &b in Nucleotide::NUCLEOTIDES.iter() {
            for &c in Nucleotide::NUCLEOTIDES.iter() {
                if a == Nucleotide::N || b == Nucleotide::N || c == Nucleotide::N {
                    let codon = Codon([a, b, c]);
                    let codon_idx: usize = CodonIdx::from(codon).into();

                    for table_idx in 0..TranslationTable::N_TRANS_TABLES {
                        let protein = ambiguous_codon_protein(codon, table_idx, translation_tables);
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
}
