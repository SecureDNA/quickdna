// Copyright 2021-2024 SecureDNA Stiftung (SecureDNA Foundation) <licensing@securedna.org>
// SPDX-License-Identifier: MIT OR Apache-2.0

use smallvec::SmallVec;

use crate::{Nucleotide, NucleotideAmbiguous, NucleotideLike};

/// Helper trait to support iters regardless of whether their items are by-ref or by-value
pub trait ToNucleotideLike
where
    // It's much easier for the compiler to reason about chained NucleotideIters
    // if we explicitly state that ToNucleotideLike is idempotent.
    Self::NucleotideType: ToNucleotideLike<NucleotideType = Self::NucleotideType>,
{
    type NucleotideType: NucleotideLike + ToNucleotideLike;

    fn to_nucleotide_like(self) -> Self::NucleotideType;
}

impl<N: NucleotideLike> ToNucleotideLike for N {
    type NucleotideType = N;

    fn to_nucleotide_like(self) -> N {
        self
    }
}

// Sadly, without specialization we can't get this to work with arbitrary
// &T where T: NucleotideLike but that shouldn't be necessary in practice.

impl ToNucleotideLike for &Nucleotide {
    type NucleotideType = Nucleotide;

    fn to_nucleotide_like(self) -> Nucleotide {
        *self
    }
}

impl ToNucleotideLike for &NucleotideAmbiguous {
    type NucleotideType = NucleotideAmbiguous;

    fn to_nucleotide_like(self) -> NucleotideAmbiguous {
        *self
    }
}

/// Extension trait for nucleotide iterators
pub trait NucleotideIter: Iterator + Sized {
    /// Returns up to 6 non-empty codon iterators for forward and reverse complement reading frames.
    ///
    /// The foward codon iterators are given before the reverse complement ones, with the forward
    /// codon iterators given in ascending order by offset from the beginning of the nucleotide
    /// sequence and the reverse complement codon iterators given in ascending order by offset
    /// from the end.
    ///
    /// # Examples
    ///
    /// ```
    /// use quickdna::{Nucleotide, NucleotideIter};
    ///
    /// use Nucleotide::*;
    /// let dna = [C, G, A, T, C, G, A, T];
    ///
    /// let frames = dna.iter().all_reading_frames();
    /// assert_eq!(frames.len(), 6);
    /// assert!(frames[0].clone().eq([
    ///     [C, G, A].into(),
    ///     [T, C, G].into(),
    /// ]));
    /// assert!(frames[1].clone().eq([
    ///     [G, A, T].into(),
    ///     [C, G, A].into(),
    /// ]));
    /// assert!(frames[2].clone().eq([
    ///     [A, T, C].into(),
    ///     [G, A, T].into(),
    /// ]));
    /// assert!(frames[3].clone().eq([
    ///     [A, T, C].into(),
    ///     [G, A, T].into(),
    /// ]));
    /// assert!(frames[4].clone().eq([
    ///     [T, C, G].into(),
    ///     [A, T, C].into(),
    /// ]));
    /// assert!(frames[5].clone().eq([
    ///     [C, G, A].into(),
    ///     [T, C, G].into(),
    /// ]));
    ///
    /// // The last forward and RC reading frames are omitted
    /// // because they would only have 2 nucleotides.
    /// let frames = dna[..4].iter().all_reading_frames();
    /// assert_eq!(frames.len(), 4);
    /// assert!(frames[0].clone().eq([
    ///     [C, G, A].into(),
    /// ]));
    /// assert!(frames[1].clone().eq([
    ///     [G, A, T].into(),
    /// ]));
    /// assert!(frames[2].clone().eq([
    ///     [A, T, C].into(),
    /// ]));
    /// assert!(frames[3].clone().eq([
    ///     [T, C, G].into(),
    /// ]));
    ///
    /// // All reading frames are omitted due to insufficient nucleotides.
    /// let frames = dna[..2].iter().all_reading_frames();
    /// assert!(frames.is_empty());
    /// ```
    fn all_reading_frames(self) -> SmallVec<[ForwardOrRcCodons<Self>; 6]>
    where
        Self: Clone + DoubleEndedIterator + ExactSizeIterator;

    /// Returns iterator of codons for the first reading frame of this nucleotide sequence.
    /// If the number of nucleotides isn't divisible by 3, excess nucleotides are silently
    /// discarded. Note that if the returned iterator is non-empty, it is the same as the
    /// first element of [`self_reading_frames`](Self::self_reading_frames).
    ///
    /// # Examples
    ///
    /// ```
    /// use quickdna::{Codons, Nucleotide, NucleotideIter};
    ///
    /// use Nucleotide::*;
    /// let dna = [C, G, A, T, C, G, A, T];
    ///
    /// let expected_codons = [
    ///     [C, G, A].into(),
    ///     [T, C, G].into(),
    /// ];
    /// assert!(dna.iter().codons().eq(expected_codons));
    ///
    /// assert!(dna.iter().codons().eq(dna.iter().self_reading_frames().remove(0)));
    /// ```
    fn codons(self) -> Codons<Self>;

    /// Returns iterator of complementary nucleotides.
    ///
    /// # Examples
    ///
    /// ```
    /// use quickdna::{Nucleotide, NucleotideIter};
    ///
    /// use Nucleotide::*;
    /// let dna = [C, G, A, T];
    ///
    /// assert!(dna.iter().complement().eq([G, C, T, A]));
    /// ```
    fn complement(self) -> Complement<Self>;

    /// Returns iterator of reverse complement of contained nucleotides.
    ///
    /// # Examples
    ///
    /// ```
    /// use quickdna::{Nucleotide, NucleotideIter};
    ///
    /// use Nucleotide::*;
    /// let dna = [C, G, A, T];
    ///
    /// assert!(dna.iter().reverse_complement().eq([A, T, C, G]));
    /// ```
    fn reverse_complement(self) -> Complement<std::iter::Rev<Self>>
    where
        Self: DoubleEndedIterator;

    /// Returns up to 3 non-empty codon iterators for reading frames.
    ///
    /// The iterators are given in ascending order of offset from the beginning of the
    /// nucleotide sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use quickdna::{Nucleotide, NucleotideIter};
    ///
    /// use Nucleotide::*;
    /// let dna = [C, G, A, T, C, G, A, T];
    ///
    /// let frames = dna.iter().self_reading_frames();
    /// assert_eq!(frames.len(), 3);
    /// assert!(frames[0].clone().eq([
    ///     [C, G, A].into(),
    ///     [T, C, G].into(),
    /// ]));
    /// assert!(frames[1].clone().eq([
    ///     [G, A, T].into(),
    ///     [C, G, A].into(),
    /// ]));
    /// assert!(frames[2].clone().eq([
    ///     [A, T, C].into(),
    ///     [G, A, T].into(),
    /// ]));
    ///
    /// // The last reading frame is omitted because it would only have 2 nucleotides.
    /// let frames = dna[..4].iter().self_reading_frames();
    /// assert_eq!(frames.len(), 2);
    /// assert!(frames[0].clone().eq([
    ///     [C, G, A].into(),
    /// ]));
    /// assert!(frames[1].clone().eq([
    ///     [G, A, T].into(),
    /// ]));
    ///
    /// // All reading frames are omitted due to insufficient nucleotides.
    /// let frames = dna[..2].iter().self_reading_frames();
    /// assert!(frames.is_empty());
    /// ```
    fn self_reading_frames(self) -> SmallVec<[Codons<Self>; 3]>
    where
        Self: Clone + ExactSizeIterator;

    /// Trims excess nucleotides off iterator end so it aligns with a codon boundary.
    ///
    /// This makes the iterator's length a multiple of 3 by removing up to 2 elements from its end.
    ///
    /// # Examples
    ///
    /// ```
    /// use quickdna::{Nucleotide, NucleotideIter};
    ///
    /// use Nucleotide::*;
    /// let dna = [C, G, A, T, C, G, A, T];
    ///
    /// let mut iter = dna.into_iter();
    /// iter.trim_to_codon();
    ///
    /// assert!(iter.eq([C, G, A, T, C, G]));
    /// ```
    fn trim_to_codon(&mut self)
    where
        Self: DoubleEndedIterator + ExactSizeIterator;

    /// Trims excess nucleotides off iterator end so it aligns with a codon boundary.
    ///
    /// This makes the iterator's length a multiple of 3 by removing up to 2 elements from its end.
    ///
    /// # Examples
    ///
    /// ```
    /// use quickdna::{Nucleotide, NucleotideIter};
    ///
    /// use Nucleotide::*;
    /// let dna = [C, G, A, T, C, G, A, T];
    ///
    /// assert!(dna.into_iter().trimmed_to_codon().eq([C, G, A, T, C, G]));
    /// ```
    fn trimmed_to_codon(mut self) -> Self
    where
        Self: DoubleEndedIterator + ExactSizeIterator,
    {
        self.trim_to_codon();
        self
    }
}

impl<N, I> NucleotideIter for I
where
    N: ToNucleotideLike,
    I: Iterator<Item = N>,
{
    fn all_reading_frames(self) -> SmallVec<[ForwardOrRcCodons<Self>; 6]>
    where
        Self: Clone + DoubleEndedIterator + ExactSizeIterator,
    {
        let iter1 = self;
        let mut iter2 = iter1.clone();
        iter2.next();
        let mut iter3 = iter2.clone();
        iter3.next();

        let iter1_rc = iter1.clone().reverse_complement();
        let mut iter2_rc = iter1_rc.clone();
        iter2_rc.next();
        let mut iter3_rc = iter2_rc.clone();
        iter3_rc.next();

        let mut frames = SmallVec::from([
            ForwardOrRcCodons::Forward(iter1.codons()),
            ForwardOrRcCodons::Forward(iter2.codons()),
            ForwardOrRcCodons::Forward(iter3.codons()),
            ForwardOrRcCodons::Rc(iter1_rc.codons()),
            ForwardOrRcCodons::Rc(iter2_rc.codons()),
            ForwardOrRcCodons::Rc(iter3_rc.codons()),
        ]);
        frames.retain(|frame| frame.len() > 0);
        frames
    }

    fn codons(self) -> Codons<Self> {
        Codons(self)
    }

    fn complement(self) -> Complement<Self> {
        Complement(self)
    }

    fn reverse_complement(self) -> Complement<std::iter::Rev<Self>>
    where
        Self: DoubleEndedIterator,
    {
        self.rev().complement()
    }

    fn self_reading_frames(self) -> SmallVec<[Codons<Self>; 3]>
    where
        Self: Clone + ExactSizeIterator,
    {
        let iter1 = self;
        let mut iter2 = iter1.clone();
        iter2.next();
        let mut iter3 = iter2.clone();
        iter3.next();
        let mut frames = SmallVec::from([iter1.codons(), iter2.codons(), iter3.codons()]);
        frames.retain(|frame| frame.len() > 0);
        frames
    }

    fn trim_to_codon(&mut self)
    where
        Self: DoubleEndedIterator + ExactSizeIterator,
    {
        for _ in 0..(self.len() % 3) {
            self.next_back();
        }
    }
}

/// Adapter yielding codons of the contained iterator.
///
/// This `struct` is created by the [`codons`](NucleotideIter::codons)
/// method on [`NucleotideIter`]. See its documentation for more.
#[derive(Clone, Debug)]
pub struct Codons<I>(I);

impl<N, I> Iterator for Codons<I>
where
    N: ToNucleotideLike,
    I: Iterator<Item = N>,
{
    type Item = <N::NucleotideType as NucleotideLike>::Codon;

    fn next(&mut self) -> Option<Self::Item> {
        match (self.0.next(), self.0.next(), self.0.next()) {
            (Some(n1), Some(n2), Some(n3)) => {
                Some([n1, n2, n3].map(|n| n.to_nucleotide_like()).into())
            }
            _ => None,
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let (min, max) = self.0.size_hint();
        (min / 3, max.map(|m| m / 3))
    }
}

impl<N, I> DoubleEndedIterator for Codons<I>
where
    N: ToNucleotideLike,
    I: DoubleEndedIterator<Item = N> + ExactSizeIterator,
{
    fn next_back(&mut self) -> Option<Self::Item> {
        self.0.trim_to_codon();
        match (self.0.next_back(), self.0.next_back(), self.0.next_back()) {
            (Some(n3), Some(n2), Some(n1)) => {
                Some([n1, n2, n3].map(|n| n.to_nucleotide_like()).into())
            }
            _ => None,
        }
    }
}

impl<I> ExactSizeIterator for Codons<I>
where
    Self: Iterator,
    I: ExactSizeIterator,
{
    fn len(&self) -> usize {
        self.0.len() / 3
    }
}

/// Adapter yielding complementary nucleotide of the contained iterator.
///
/// This `struct` is created by the [`complement`](NucleotideIter::complement)
/// method on [`NucleotideIter`]. See its documentation for more.
#[derive(Clone, Debug)]
pub struct Complement<I>(I);

impl<N, I> Iterator for Complement<I>
where
    N: ToNucleotideLike,
    I: Iterator<Item = N>,
{
    type Item = N::NucleotideType;

    fn next(&mut self) -> Option<Self::Item> {
        self.0.next().map(|n| n.to_nucleotide_like().complement())
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.0.size_hint()
    }
}

impl<N, I> DoubleEndedIterator for Complement<I>
where
    N: ToNucleotideLike,
    I: DoubleEndedIterator<Item = N>,
{
    fn next_back(&mut self) -> Option<Self::Item> {
        self.0
            .next_back()
            .map(|n| n.to_nucleotide_like().complement())
    }
}

impl<I> ExactSizeIterator for Complement<I>
where
    Self: Iterator,
    I: ExactSizeIterator,
{
    fn len(&self) -> usize {
        self.0.len()
    }
}

/// Adapter capable of holding either forward codon iterators or reverse complement codon iterators.
///
/// This `struct` is created by the [`all_reading_frames`](NucleotideIter::all_reading_frames)
/// method on [`NucleotideIter`]. See its documentation for more.
#[derive(Clone, Debug)]
pub enum ForwardOrRcCodons<I> {
    Forward(Codons<I>),
    Rc(Codons<Complement<std::iter::Rev<I>>>),
}

impl<N, I> Iterator for ForwardOrRcCodons<I>
where
    N: ToNucleotideLike,
    I: DoubleEndedIterator<Item = N>,
{
    type Item = <N::NucleotideType as NucleotideLike>::Codon;

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            Self::Forward(iter) => iter.next(),
            Self::Rc(iter) => iter.next(),
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        match self {
            Self::Forward(iter) => iter.size_hint(),
            Self::Rc(iter) => iter.size_hint(),
        }
    }
}

impl<N, I> DoubleEndedIterator for ForwardOrRcCodons<I>
where
    N: ToNucleotideLike,
    I: DoubleEndedIterator<Item = N> + ExactSizeIterator,
{
    fn next_back(&mut self) -> Option<Self::Item> {
        match self {
            Self::Forward(iter) => iter.next_back(),
            Self::Rc(iter) => iter.next_back(),
        }
    }
}

impl<N, I> ExactSizeIterator for ForwardOrRcCodons<I>
where
    Self: Iterator,
    N: ToNucleotideLike,
    I: DoubleEndedIterator<Item = N> + ExactSizeIterator,
{
    fn len(&self) -> usize {
        match self {
            Self::Forward(iter) => iter.len(),
            Self::Rc(iter) => iter.len(),
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_reverse_codons() {
        use Nucleotide::*;
        let dna = [A, A, T, T, C, C, G, G];
        let rev_codons: Vec<_> = dna.iter().codons().rev().collect();
        let expected = [[T, C, C].into(), [A, A, T].into()];
        assert_eq!(rev_codons, expected);
    }

    #[test]
    fn can_compile_iter_of_generic_nucleotide_references() {
        // The previous version of `NucleotideIter` worked with iterators of `Nucleotide`,
        // `&Nucleotide`, `NucleotideAmbiguous`, `&NucleotideAmbiguous` and `impl NucleotideLike`,
        // but NOT `&impl NucleotideLike`.
        fn _do_stuff(dna: &[impl NucleotideLike]) {
            dna.iter().reverse_complement(); // no-op; just checking that it compiles
        }
    }
}
