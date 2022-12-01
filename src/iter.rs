use smallvec::SmallVec;

use crate::{Nucleotide, NucleotideAmbiguous, NucleotideLike};

/// Helper trait to support iters regardless of whether their items are by-ref or by-value
pub trait ToNucleotideLike
where
    // It's much easier for the compiler to reason about chained Nucleotides iterator
    // adapters if we explicitly state that ToNucleotideLike is idempotent.
    Self::NucleotideType: ToNucleotideLike<NucleotideType = Self::NucleotideType>,
{
    type NucleotideType: NucleotideLike + ToNucleotideLike;

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
    /// use quickdna::{Nucleotide, Nucleotides};
    ///
    /// use Nucleotide::*;
    /// let dna = [C, G, A, T, C, G, A, T];
    ///
    /// let frames = dna.all_reading_frames();
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
    /// let frames = dna[..4].all_reading_frames();
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
    /// let frames = dna[..2].all_reading_frames();
    /// assert!(frames.is_empty());
    /// ```
    fn all_reading_frames(self) -> SmallVec<[ForwardOrRcCodons<Self::IntoIter>; 6]>
    where
        Self::IntoIter: Clone + DoubleEndedIterator + ExactSizeIterator;

    /// Returns iterator of codons for the first reading frame of this nucleotide sequence.
    /// If the number of nucleotides isn't divisible by 3, excess nucleotides are silently
    /// discarded. Note that if the returned iterator is non-empty, it is the same as the
    /// first element of [`self_reading_frames`](Self::self_reading_frames).
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
    /// assert!(dna.codons().eq(dna.self_reading_frames().remove(0)));
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
    /// let frames = dna.self_reading_frames();
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
    /// let frames = dna[..4].self_reading_frames();
    /// assert_eq!(frames.len(), 2);
    /// assert!(frames[0].clone().eq([
    ///     [C, G, A].into(),
    /// ]));
    /// assert!(frames[1].clone().eq([
    ///     [G, A, T].into(),
    /// ]));
    ///
    /// // All reading frames are omitted due to insufficient nucleotides.
    /// let frames = dna[..2].self_reading_frames();
    /// assert!(frames.is_empty());
    /// ```
    fn self_reading_frames(self) -> SmallVec<[Codons<Self::IntoIter>; 3]>
    where
        Self::IntoIter: Clone + ExactSizeIterator;
}

impl<N, I, T> Nucleotides for T
where
    N: ToNucleotideLike,
    I: Iterator<Item = N>,
    T: IntoIterator<IntoIter = I>,
{
    fn all_reading_frames(self) -> SmallVec<[ForwardOrRcCodons<Self::IntoIter>; 6]>
    where
        Self::IntoIter: Clone + DoubleEndedIterator + ExactSizeIterator,
    {
        let iter1 = self.into_iter();
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

    fn codons(self) -> Codons<Self::IntoIter> {
        Codons(self.into_iter())
    }

    fn complement(self) -> Complement<Self::IntoIter> {
        Complement(self.into_iter())
    }

    fn reverse_complement(self) -> std::iter::Rev<Complement<Self::IntoIter>>
    where
        Self::IntoIter: DoubleEndedIterator,
    {
        self.complement().rev()
    }

    fn self_reading_frames(self) -> SmallVec<[Codons<Self::IntoIter>; 3]>
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
    I: DoubleEndedIterator<Item = N> + ExactSizeIterator,
{
    fn next_back(&mut self) -> Option<Self::Item> {
        let remainder_nucleotides = self.0.len() % 3;
        for _ in 0..remainder_nucleotides {
            self.0.next_back();
        }
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
    fn len(&self) -> usize {
        self.0.len()
    }
}

/// Adapter capable of holding either forward codon iterators or reverse complement codon iterators.
///
/// This `struct` is created by the [`all_reading_frames`](Nucleotides::all_reading_frames)
/// method on [`Nucleotides`]. See its documentation for more.
#[derive(Clone, Debug)]
pub enum ForwardOrRcCodons<I> {
    Forward(Codons<I>),
    Rc(Codons<std::iter::Rev<Complement<I>>>),
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
        let rev_codons: Vec<_> = dna.codons().rev().collect();
        let expected = [[T, C, C].into(), [A, A, T].into()];
        assert_eq!(rev_codons, expected);
    }
}
