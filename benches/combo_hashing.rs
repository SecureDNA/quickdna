// Copyright 2021-2023 SecureDNA Stiftung (SecureDNA Foundation) <licensing@securedna.org>
// SPDX-License-Identifier: MIT OR Apache-2.0

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use rand::{rngs::OsRng, seq::SliceRandom};
// use smallvec::SmallVec;

use quickdna::permute::{canonicalize, LexicalMin};
use quickdna::Nucleotide;

struct ForwardCanonicalized2<I> {
    inner: I,
    permutation: [Option<Nucleotide>; 4],
    unmapped: std::slice::Iter<'static, Nucleotide>,
}

impl<I> ForwardCanonicalized2<I>
where
    I: Iterator<Item = Nucleotide>,
{
    fn new(iterable: impl IntoIterator<IntoIter = I>) -> Self {
        Self {
            inner: iterable.into_iter(),
            permutation: [None; 4],
            unmapped: Nucleotide::ALL.iter(),
        }
    }
}

impl<I> Iterator for ForwardCanonicalized2<I>
where
    I: Iterator<Item = Nucleotide>,
{
    type Item = Nucleotide;

    fn next(&mut self) -> Option<Self::Item> {
        self.inner.next().map(|nuc| {
            let idx = match nuc {
                Nucleotide::A => 0,
                Nucleotide::T => 1,
                Nucleotide::C => 2,
                Nucleotide::G => 3,
            };
            *self.permutation[idx]
                .get_or_insert_with(|| self.unmapped.next().copied().unwrap_or(Nucleotide::A))
        })
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.inner.size_hint()
    }
}

impl<I: ExactSizeIterator<Item = Nucleotide>> ExactSizeIterator for ForwardCanonicalized2<I> {}

pub struct Canonicalized2<I>(
    LexicalMin<ForwardCanonicalized2<I>, ForwardCanonicalized2<std::iter::Rev<I>>>,
);

impl<I> Iterator for Canonicalized2<I>
where
    I: DoubleEndedIterator<Item = Nucleotide>,
{
    type Item = Nucleotide;

    fn next(&mut self) -> Option<Self::Item> {
        self.0.next()
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.0.size_hint()
    }
}

impl<I> ExactSizeIterator for Canonicalized2<I> where
    I: ExactSizeIterator + DoubleEndedIterator<Item = Nucleotide>
{
}

fn canonicalize2<I>(iterable: impl IntoIterator<IntoIter = I>) -> Canonicalized2<I>
where
    I: DoubleEndedIterator<Item = Nucleotide> + Clone,
{
    let iter = iterable.into_iter();
    let rev_canon = ForwardCanonicalized2::new(iter.clone().rev());
    let fw_canon = ForwardCanonicalized2::new(iter);
    Canonicalized2(LexicalMin::new(fw_canon, rev_canon))
}

// use quickdna::NucleotideLike;
//
// struct ForwardCanonicalized2<I> {
//     inner: I,
//     is_mapped: u8,
//     permutation: [Nucleotide; 4],
//     unmapped: u8,
// }
//
// impl<I> ForwardCanonicalized2<I>
// where
//     I: Iterator<Item = Nucleotide>,
// {
//     fn new(iterable: impl IntoIterator<IntoIter = I>) -> Self {
//         Self {
//             inner: iterable.into_iter(),
//             is_mapped: 0,
//             permutation: [Nucleotide::A; 4],
//             unmapped: 3,
//         }
//     }
// }
//
// impl<I> Iterator for ForwardCanonicalized2<I>
// where
//     I: Iterator<Item = Nucleotide>,
// {
//     type Item = Nucleotide;
//
//     fn next(&mut self) -> Option<Self::Item> {
//         self.inner.next().map(|nuc| {
//             let bits = nuc.bits();
//             let mapped_nuc = &mut self.permutation[bits.ilog2() as usize];
//             if bits & self.is_mapped == 0 {
//                 self.is_mapped |= bits;
//                 self.unmapped = (self.unmapped + 1) % 4;
//                 *mapped_nuc = Nucleotide::ALL[self.unmapped as usize];
//             }
//             *mapped_nuc
//         })
//     }
//
//     fn size_hint(&self) -> (usize, Option<usize>) {
//         self.inner.size_hint()
//     }
// }
//
// impl<I: ExactSizeIterator<Item = Nucleotide>> ExactSizeIterator for ForwardCanonicalized2<I> {}
//
// pub struct Canonicalized2<I>(
//     LexicalMin<ForwardCanonicalized2<I>, ForwardCanonicalized2<std::iter::Rev<I>>>,
// );
//
// impl<I> Iterator for Canonicalized2<I>
// where
//     I: DoubleEndedIterator<Item = Nucleotide>,
// {
//     type Item = Nucleotide;
//
//     fn next(&mut self) -> Option<Self::Item> {
//         self.0.next()
//     }
//
//     fn size_hint(&self) -> (usize, Option<usize>) {
//         self.0.size_hint()
//     }
// }
//
// impl<I> ExactSizeIterator for Canonicalized2<I> where
//     I: ExactSizeIterator + DoubleEndedIterator<Item = Nucleotide>
// {
// }
//
// fn canonicalize2<I>(iterable: impl IntoIterator<IntoIter = I>) -> Canonicalized2<I>
// where
//     I: DoubleEndedIterator<Item = Nucleotide> + Clone,
// {
//     let iter = iterable.into_iter();
//     let rev_canon = ForwardCanonicalized2::new(iter.clone().rev());
//     let fw_canon = ForwardCanonicalized2::new(iter);
//     Canonicalized2(LexicalMin::new(fw_canon, rev_canon))
// }

pub fn criterion_benchmark(c: &mut Criterion) {
    const NUM_WINDOWS: u64 = 10000;
    const WINDOW_LEN: usize = 42;

    let mut group = c.benchmark_group("2-way canonicalization");
    group.throughput(Throughput::Elements(NUM_WINDOWS));
    group.bench_with_input(
        BenchmarkId::from_parameter(NUM_WINDOWS),
        &NUM_WINDOWS,
        |b, &num_windows| {
            let dna: Vec<_> = (0..num_windows as usize + WINDOW_LEN - 1)
                .map(|_| Nucleotide::ALL.choose(&mut OsRng).unwrap())
                .collect();

            b.iter(|| {
                let canonicalized_windows = dna.windows(WINDOW_LEN).map(|window| {
                    let mut canonical = canonicalize(window.iter().copied().copied());
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
    drop(group);

    let mut group = c.benchmark_group("2-way canonicalization unambiguous");
    group.throughput(Throughput::Elements(NUM_WINDOWS));
    group.bench_with_input(
        BenchmarkId::from_parameter(NUM_WINDOWS),
        &NUM_WINDOWS,
        |b, &num_windows| {
            let dna: Vec<_> = (0..num_windows as usize + WINDOW_LEN - 1)
                .map(|_| Nucleotide::ALL.choose(&mut OsRng).unwrap())
                .collect();

            b.iter(|| {
                let canonicalized_windows = dna.windows(WINDOW_LEN).map(|window| {
                    let mut canonical = canonicalize2(window.iter().copied().copied());
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
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
