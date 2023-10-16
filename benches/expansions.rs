// Copyright 2021-2023 SecureDNA Stiftung (SecureDNA Foundation) <licensing@securedna.org>
// SPDX-License-Identifier: MIT OR Apache-2.0

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use quickdna::{
    expansions::Expansions, BaseSequence, DnaSequenceStrict, Nucleotide, NucleotideAmbiguous,
};
use rand::{rngs::OsRng, seq::SliceRandom};

// Vec-based expansions code taken from DNA screening tools

const MAX_EXPANSIONS_FOR_WINDOW: usize = 16;

fn vec_based_expansions<const WINDOW_SIZE: usize>(
    wildcards: [NucleotideAmbiguous; WINDOW_SIZE],
) -> Vec<[Nucleotide; WINDOW_SIZE]> {
    if !are_expansions_feasible(wildcards) {
        panic!(
            "Too many possible expansions of window {:?}; ignoring window",
            wildcards
        );
    }

    let mut expanded_sequence_buffers: Vec<DnaSequenceStrict> =
        vec![DnaSequenceStrict::new(vec![])];

    for wildcard in wildcards {
        let expansions = wildcard.possibilities();
        if expansions.len() == 1 {
            for sequence_buffer in expanded_sequence_buffers.iter_mut() {
                sequence_buffer.push(expansions[0])
            }
        } else {
            let mut new_buffers = Vec::new();
            for nucleotide_from_expansion in expansions.iter() {
                for existing_sequence_buffer in expanded_sequence_buffers.iter() {
                    let mut new_sequence_bufer = existing_sequence_buffer.clone();
                    new_sequence_bufer.push(*nucleotide_from_expansion);
                    new_buffers.push(new_sequence_bufer);
                }
            }
            expanded_sequence_buffers = new_buffers;
        }
    }

    expanded_sequence_buffers
        .into_iter()
        .map(|s| s.as_slice().try_into().unwrap())
        .collect()
}

/// each expansion creates n variants
/// for example AAR becomes [AAA, AAG]
/// these expansions are exponential
/// AARR becomes [AAAA, AAAG, AAGA, AAGG] and so on
/// we do not want to generate more than a limited number of these
fn are_expansions_feasible(ws: impl IntoIterator<Item = NucleotideAmbiguous>) -> bool {
    let mut acc = 1;
    for w in ws {
        acc *= w.possibilities().len();
        if acc > MAX_EXPANSIONS_FOR_WINDOW {
            return false;
        }
    }
    acc <= MAX_EXPANSIONS_FOR_WINDOW
}

fn semi_ambiguous_dna(dna_len: usize, num_ambiguities: usize) -> Vec<NucleotideAmbiguous> {
    let nucleotides = Nucleotide::ALL.map(|nuc| nuc.into());
    let ambiguous_nucleotides = [
        NucleotideAmbiguous::W,
        NucleotideAmbiguous::M,
        NucleotideAmbiguous::Y,
        NucleotideAmbiguous::H,
        NucleotideAmbiguous::R,
        NucleotideAmbiguous::K,
        NucleotideAmbiguous::D,
        NucleotideAmbiguous::S,
        NucleotideAmbiguous::V,
        NucleotideAmbiguous::B,
        NucleotideAmbiguous::N,
    ];
    let mut dna: Vec<_> = (0..dna_len)
        .map(|_| *nucleotides.choose(&mut OsRng).unwrap())
        .collect();
    for i in rand::seq::index::sample(&mut OsRng, dna_len, num_ambiguities) {
        dna[i] = *ambiguous_nucleotides.choose(&mut OsRng).unwrap();
    }
    dna
}

pub fn criterion_benchmark(c: &mut Criterion) {
    let num_ambiguities = 2;
    const WINDOW_LEN: usize = 42;
    let windows: Vec<[_; WINDOW_LEN]> = (0..1000)
        .map(|_| {
            semi_ambiguous_dna(WINDOW_LEN, num_ambiguities)
                .try_into()
                .unwrap()
        })
        .collect();

    c.bench_function("Vec-based ambiguity expansions", |b| {
        b.iter(|| {
            for window in &windows {
                for expansion in vec_based_expansions(*window) {
                    black_box(expansion);
                }
            }
        })
    });

    c.bench_function("Iter-based ambiguity expansions", |b| {
        b.iter(|| {
            for window in &windows {
                for expansion in Expansions::new(window) {
                    black_box(expansion);
                }
            }
        })
    });

    c.bench_function("Iter-based number of expansions", |b| {
        b.iter(|| {
            for window in &windows {
                black_box(Expansions::new(window).size_hint());
            }
        })
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
