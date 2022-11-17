use std::iter;
use std::slice;

use smallvec::SmallVec;

pub use crate::amino_acid::AminoAcid;
use crate::NucleotideLike;

pub trait Nucleotides<N: NucleotideLike>
where
    Self: Sized,
{
    type Codons;
    type ReverseComplement;

    fn codons(self) -> Self::Codons;
    fn frames(self) -> SmallVec<[Self; 3]>
    where
        Self: Clone;
    fn reverse_complement(self) -> Self::ReverseComplement;
}

impl<'a, N: NucleotideLike> Nucleotides<N> for &'a [N] {
    type Codons = Codons<iter::Copied<slice::Iter<'a, N>>>;
    type ReverseComplement = ReverseComplement<iter::Copied<slice::Iter<'a, N>>>;

    fn codons(self) -> Self::Codons {
        self.iter().copied().codons()
    }

    fn frames(self) -> SmallVec<[Self; 3]> {
        let mut v = SmallVec::from_buf([&self[0..], &self[1..], &self[2..]]);
        v.retain(|f| f.len() >= 3);
        v
    }

    fn reverse_complement(self) -> Self::ReverseComplement {
        self.iter().copied().reverse_complement()
    }
}

impl<'a, N: NucleotideLike> Nucleotides<N> for iter::Copied<slice::Iter<'a, N>> {
    type Codons = Codons<Self>;
    type ReverseComplement = ReverseComplement<Self>;

    fn codons(self) -> Self::Codons {
        Codons(self)
    }

    fn frames(self) -> SmallVec<[Self; 3]> {
        frames(self)
    }

    fn reverse_complement(self) -> Self::ReverseComplement {
        ReverseComplement(self)
    }
}

impl<N: NucleotideLike, I: DoubleEndedIterator<Item = N> + ExactSizeIterator> Nucleotides<N>
    for ReverseComplement<I>
{
    type Codons = Codons<Self>;
    type ReverseComplement = I;

    fn codons(self) -> Self::Codons {
        Codons(self)
    }

    fn frames(self) -> SmallVec<[Self; 3]>
    where
        Self: Clone,
    {
        frames(self)
    }

    fn reverse_complement(self) -> Self::ReverseComplement {
        self.0
    }
}

fn frames<I: ExactSizeIterator + Clone>(mut iter: I) -> SmallVec<[I; 3]> {
    let frame0 = iter.clone();
    iter.next();
    let frame1 = iter.clone();
    iter.next();
    let frame2 = iter;
    let mut frames = SmallVec::from_buf([frame0, frame1, frame2]);
    frames.retain(|f| f.len() >= 3);
    frames
}

// We create this type explicitly (as opposed to via .rev().map(...))
// so that it has a concrete name, and so that it's easy to peel off
#[derive(Clone)]
pub struct ReverseComplement<I: DoubleEndedIterator>(I);

impl<N: NucleotideLike, I: DoubleEndedIterator<Item = N>> Iterator for ReverseComplement<I> {
    type Item = N;

    fn next(&mut self) -> Option<N> {
        self.0.next_back().map(|n| n.complement())
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.0.size_hint()
    }
}

impl<N: NucleotideLike, I: DoubleEndedIterator<Item = N>> DoubleEndedIterator
    for ReverseComplement<I>
{
    fn next_back(&mut self) -> Option<N> {
        self.0.next().map(|n| n.complement())
    }
}

impl<N: NucleotideLike, I: DoubleEndedIterator<Item = N> + ExactSizeIterator> ExactSizeIterator
    for ReverseComplement<I>
{
}

#[derive(Clone)]
pub struct Codons<I: Iterator>(I);

impl<N: NucleotideLike, I: Iterator<Item = N>> Iterator for Codons<I> {
    type Item = N::Codon;

    fn next(&mut self) -> Option<N::Codon> {
        match (self.0.next(), self.0.next(), self.0.next()) {
            (Some(n1), Some(n2), Some(n3)) => Some([n1, n2, n3].into()),
            _ => None,
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let (min, max) = self.0.size_hint();
        (min / 3, max.map(|m| m / 3))
    }
}

impl<N: NucleotideLike, I: ExactSizeIterator<Item = N>> ExactSizeIterator for Codons<I> {}

impl<N: NucleotideLike, I: DoubleEndedIterator<Item = N> + ExactSizeIterator> DoubleEndedIterator
    for Codons<I>
{
    fn next_back(&mut self) -> Option<Self::Item> {
        for _ in 0..self.0.len() % 3 {
            let _ = self.0.next_back();
        }
        match (self.0.next_back(), self.0.next_back(), self.0.next_back()) {
            (Some(n1), Some(n2), Some(n3)) => Some([n3, n2, n1].into()),
            _ => None,
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    use std::str::FromStr;

    use crate::{Codon, DnaSequence, Nucleotide, NucleotideAmbiguous, TranslationTable};

    #[test]
    fn example_code() {
        use Nucleotide::*;

        let dna: DnaSequence<NucleotideAmbiguous> = DnaSequence::from_str("ATCGAATTCCGG").unwrap();
        let codons: Vec<_> = dna.codons().collect();
        let expected = [
            Codon::from([A, T, C]).into(),
            Codon::from([G, A, A]).into(),
            Codon::from([T, T, C]).into(),
            Codon::from([C, G, G]).into(),
        ];
        assert_eq!(codons, expected);

        let array = [A, T, C, G, A, A, T, T, C, C, G, G];
        let rc_frame1: Vec<_> = (&array).reverse_complement().frames()[1]
            .clone()
            .codons()
            .collect();
        let expected = [
            Codon::from([C, G, G]),
            Codon::from([A, A, T]),
            Codon::from([T, C, G]),
        ];
        assert_eq!(rc_frame1, expected);

        let slice = &array[3..];
        let maybe_protein: Option<SmallVec<[AminoAcid; 10]>> = slice
            .codons()
            .map(TranslationTable::Ncbi1.to_fn())
            .collect();
        let protein = maybe_protein.unwrap();
        assert_eq!(
            protein.as_slice(),
            &[AminoAcid::E, AminoAcid::F, AminoAcid::R]
        );

        // Probably not useful, but this should be about as fast as iterating over the original array:
        let unchanged: Vec<_> = array
            .reverse_complement()
            .reverse_complement()
            .reverse_complement()
            .reverse_complement()
            .reverse_complement()
            .reverse_complement()
            .reverse_complement()
            .reverse_complement()
            .collect();
        assert_eq!(unchanged, array);
    }
}
