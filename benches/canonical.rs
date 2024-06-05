// Copyright 2021-2024 SecureDNA Stiftung (SecureDNA Foundation) <licensing@securedna.org>
// SPDX-License-Identifier: MIT OR Apache-2.0

use std::hint::black_box;

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use rand::{rngs::OsRng, seq::SliceRandom};

use quickdna::canonical::{Canonical, ForwardCanonical};
use quickdna::Nucleotide;

pub fn criterion_benchmark(c: &mut Criterion) {
    const NUM_WINDOWS: u64 = 10000;
    const WINDOW_LEN: usize = 42;
    let dna_len = (NUM_WINDOWS as usize) + WINDOW_LEN - 1;
    let dna: Vec<_> = (0..dna_len)
        .map(|_| Nucleotide::ALL.choose(&mut OsRng).unwrap())
        .collect();

    let num_windows_desc = format!("{NUM_WINDOWS} windows");

    let mut group = c.benchmark_group("canonicalization");
    group.throughput(Throughput::Elements(NUM_WINDOWS));
    group.bench_with_input(
        BenchmarkId::new("bidirectional", &num_windows_desc),
        &(&dna, WINDOW_LEN),
        |b, &(dna, window_len)| {
            b.iter(|| {
                let canonicalized_windows = dna.windows(window_len).map(|window| {
                    let mut canonical = Canonical::new(window.iter().copied().copied());
                    let window: [_; WINDOW_LEN] =
                        std::array::from_fn(|_| canonical.next().unwrap());
                    window
                });
                for window in canonicalized_windows {
                    black_box(&window);
                }
            })
        },
    );
    group.bench_with_input(
        BenchmarkId::new("forward-only", &num_windows_desc),
        &(&dna, WINDOW_LEN),
        |b, &(dna, window_len)| {
            b.iter(|| {
                let canonicalized_windows = dna.windows(window_len).map(|window| {
                    let mut canonical = ForwardCanonical::new(window.iter().copied().copied());
                    let window: [_; WINDOW_LEN] =
                        std::array::from_fn(|_| canonical.next().unwrap());
                    window
                });
                for window in canonicalized_windows {
                    black_box(&window);
                }
            })
        },
    );
    group.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
