use smallvec::SmallVec;

use crate::{Nucleotide, NucleotideAmbiguous, NucleotideLike};

/// Helper trait to support iters regardless of whether their items are by-ref or by-value
pub trait ToNucleotideLike {
    type NucleotideType: NucleotideLike;

    fn to_nucleotide_like(self) -> Self::NucleotideType;
}

impl ToNucleotideLike for Nucleotide {
    type NucleotideType = Nucleotide;

    fn to_nucleotide_like(self) -> Nucleotide {
        self
    }
}

impl ToNucleotideLike for &Nucleotide {
    type NucleotideType = Nucleotide;

    fn to_nucleotide_like(self) -> Nucleotide {
        *self
    }
}

impl ToNucleotideLike for NucleotideAmbiguous {
    type NucleotideType = NucleotideAmbiguous;

    fn to_nucleotide_like(self) -> NucleotideAmbiguous {
        self
    }
}

impl ToNucleotideLike for &NucleotideAmbiguous {
    type NucleotideType = NucleotideAmbiguous;

    fn to_nucleotide_like(self) -> NucleotideAmbiguous {
        *self
    }
}

/// Extension trait for nucleotide sequences
pub trait Nucleotides: IntoIterator {
    /// Returns iterator of codons for the first reading frame of this nucleotide sequence.
    /// If the number of nucleotides isn't divisible by 3, excess nucleotides are silently
    /// discarded. Note that if the returned iterator is non-empty, it is the same as the
    /// first element of [`reading_frames`](Self::reading_frames).
    ///
    /// # Examples
    ///
    /// ```
    /// use quickdna::{Codons, Nucleotide, Nucleotides};
    ///
    /// use Nucleotide::*;
    /// let dna = [C, G, A, T, C, G, A, T];
    ///
    /// let expected_codons = [
    ///     [C, G, A].into(),
    ///     [T, C, G].into(),
    /// ];
    /// assert!(dna.codons().eq(expected_codons));
    ///
    /// assert!(dna.codons().eq(dna.reading_frames().remove(0)));
    /// ```
    fn codons(self) -> Codons<Self::IntoIter>;

    /// Returns iterator of complementary nucleotides.
    ///
    /// # Examples
    ///
    /// ```
    /// use quickdna::{Nucleotide, Nucleotides};
    ///
    /// use Nucleotide::*;
    /// let dna = [C, G, A, T];
    ///
    /// assert!(dna.complement().eq([G, C, T, A]));
    /// ```
    fn complement(self) -> Complement<Self::IntoIter>;

    /// Returns up to 3 non-empty codon iterators for reading frames.
    ///
    /// The iterators are given in ascending order of offset from the beginning of the
    /// nucleotide sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use quickdna::{Nucleotide, Nucleotides};
    ///
    /// use Nucleotide::*;
    /// let dna = [C, G, A, T, C, G, A, T];
    ///
    /// let frames = dna.reading_frames();
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
    /// let frames = dna[..4].reading_frames();
    /// assert_eq!(frames.len(), 2);
    /// assert!(frames[0].clone().eq([
    ///     [C, G, A].into(),
    /// ]));
    /// assert!(frames[1].clone().eq([
    ///     [G, A, T].into(),
    /// ]));
    ///
    /// // All reading frames are omitted due to insufficient nucleotides.
    /// let frames = dna[..2].reading_frames();
    /// assert!(frames.is_empty());
    /// ```
    fn reading_frames(self) -> SmallVec<[Codons<Self::IntoIter>; 3]>
    where
        Self::IntoIter: Clone + ExactSizeIterator;

    /// Returns iterator of reverse complement of contained nucleotides.
    ///
    /// # Examples
    ///
    /// ```
    /// use quickdna::{Nucleotide, Nucleotides};
    ///
    /// use Nucleotide::*;
    /// let dna = [C, G, A, T];
    ///
    /// assert!(dna.reverse_complement().eq([A, T, C, G]));
    /// ```
    fn reverse_complement(self) -> std::iter::Rev<Complement<Self::IntoIter>>
    where
        Self::IntoIter: DoubleEndedIterator;
}

impl<N, I, T> Nucleotides for T
where
    N: ToNucleotideLike,
    I: Iterator<Item = N>,
    T: IntoIterator<IntoIter = I>,
{
    fn codons(self) -> Codons<Self::IntoIter> {
        Codons(self.into_iter())
    }

    fn complement(self) -> Complement<Self::IntoIter> {
        Complement(self.into_iter())
    }

    fn reading_frames(self) -> SmallVec<[Codons<Self::IntoIter>; 3]>
    where
        Self::IntoIter: Clone + ExactSizeIterator,
    {
        let iter1 = self.into_iter();
        let mut iter2 = iter1.clone();
        iter2.next();
        let mut iter3 = iter2.clone();
        iter3.next();
        let mut frames = SmallVec::from([iter1.codons(), iter2.codons(), iter3.codons()]);
        frames.retain(|frame| frame.len() > 0);
        frames
    }

    fn reverse_complement(self) -> std::iter::Rev<Complement<Self::IntoIter>>
    where
        Self::IntoIter: DoubleEndedIterator,
    {
        self.complement().rev()
    }
}

/// Adapter yielding codons of the contained iterator.
///
/// This `struct` is created by the [`codons`](Nucleotides::codons)
/// method on [`Nucleotides`]. See its documentation for more.
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
    I: DoubleEndedIterator<Item = N>,
{
    fn next_back(&mut self) -> Option<Self::Item> {
        match (self.0.next_back(), self.0.next_back(), self.0.next_back()) {
            (Some(n1), Some(n2), Some(n3)) => {
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
}

/// Adapter yielding complementary nucleotide of the contained iterator.
///
/// This `struct` is created by the [`complement`](Nucleotides::complement)
/// method on [`Nucleotides`]. See its documentation for more.
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
}
