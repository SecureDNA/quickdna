use smallvec::SmallVec;

pub use crate::amino_acid::AminoAcid;
use crate::{Nucleotide, NucleotideAmbiguous, NucleotideLike};

pub trait ToNucleotideValue {
    type NucleotideValue: NucleotideLike;

    fn to_nucleotide(self) -> Self::NucleotideValue;
}

impl ToNucleotideValue for Nucleotide {
    type NucleotideValue = Nucleotide;

    fn to_nucleotide(self) -> Nucleotide {
        self
    }
}

impl ToNucleotideValue for &Nucleotide {
    type NucleotideValue = Nucleotide;

    fn to_nucleotide(self) -> Nucleotide {
        *self
    }
}

impl ToNucleotideValue for NucleotideAmbiguous {
    type NucleotideValue = Self;

    fn to_nucleotide(self) -> NucleotideAmbiguous {
        self
    }
}

impl ToNucleotideValue for &NucleotideAmbiguous {
    type NucleotideValue = NucleotideAmbiguous;

    fn to_nucleotide(self) -> NucleotideAmbiguous {
        *self
    }
}

trait Nucleotides {
    type Codons;
    type ReverseComplement
    where
        Self: IntoIterator,
        Self::IntoIter: DoubleEndedIterator;

    fn codons(self) -> Self::Codons;

    fn frames(self) -> SmallVec<[<Self as IntoIterator>::IntoIter; 3]>
    where
        Self: IntoIterator,
        Self::IntoIter: ExactSizeIterator + Clone;

    fn reverse_complement(self) -> Self::ReverseComplement
    where
        Self: IntoIterator,
        Self::IntoIter: DoubleEndedIterator;
}

impl<I> Nucleotides for I
where
    I: IntoIterator,
    <I as IntoIterator>::Item: ToNucleotideValue,
    <<I as IntoIterator>::Item as ToNucleotideValue>::NucleotideValue: NucleotideLike,
{
    type Codons = Codons<<Self as IntoIterator>::IntoIter>;
    type ReverseComplement = ReverseComplement<<Self as IntoIterator>::IntoIter>
        where Self: IntoIterator, <Self as IntoIterator>::IntoIter: DoubleEndedIterator;

    fn codons(self) -> Self::Codons {
        Codons(self.into_iter())
    }

    fn frames(self) -> SmallVec<[<Self as IntoIterator>::IntoIter; 3]>
    where
        <Self as IntoIterator>::IntoIter: ExactSizeIterator + Clone,
    {
        let mut iter = self.into_iter();
        let frame0 = iter.clone();
        iter.next();
        let frame1 = iter.clone();
        iter.next();
        let frame2 = iter;
        let mut frames = SmallVec::from_buf([frame0, frame1, frame2]);
        frames.retain(|f| f.len() >= 3);
        frames
    }

    fn reverse_complement(self) -> Self::ReverseComplement
    where
        Self: IntoIterator,
        <Self as IntoIterator>::IntoIter: DoubleEndedIterator,
    {
        ReverseComplement(self.into_iter())
    }
}

#[derive(Clone)]
pub struct Codons<I: Iterator>(I);

impl<N: ToNucleotideValue, I: Iterator<Item = N>> Iterator for Codons<I> {
    type Item = <<N as ToNucleotideValue>::NucleotideValue as NucleotideLike>::Codon;

    fn next(&mut self) -> Option<Self::Item> {
        match (self.0.next(), self.0.next(), self.0.next()) {
            (Some(n1), Some(n2), Some(n3)) => {
                Some([n1.to_nucleotide(), n2.to_nucleotide(), n3.to_nucleotide()].into())
            }
            _ => None,
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let (min, max) = self.0.size_hint();
        (min / 3, max.map(|m| m / 3))
    }
}

impl<N: ToNucleotideValue, I: ExactSizeIterator<Item = N>> ExactSizeIterator for Codons<I> {}

impl<N: ToNucleotideValue, I: DoubleEndedIterator<Item = N> + ExactSizeIterator> DoubleEndedIterator
    for Codons<I>
{
    fn next_back(&mut self) -> Option<Self::Item> {
        for _ in 0..self.0.len() % 3 {
            let _ = self.0.next_back();
        }
        match (self.0.next_back(), self.0.next_back(), self.0.next_back()) {
            (Some(n1), Some(n2), Some(n3)) => {
                Some([n3.to_nucleotide(), n2.to_nucleotide(), n1.to_nucleotide()].into())
            }
            _ => None,
        }
    }
}

// We create this type explicitly (as opposed to via .rev().map(...))
// so that it has a concrete name, and so that it's easy to peel off
#[derive(Clone)]
pub struct ReverseComplement<I: DoubleEndedIterator>(I);

impl<N: ToNucleotideValue, I: DoubleEndedIterator<Item = N>> Iterator for ReverseComplement<I> {
    type Item = N::NucleotideValue;

    fn next(&mut self) -> Option<Self::Item> {
        self.0.next_back().map(|n| n.to_nucleotide().complement())
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.0.size_hint()
    }
}

impl<N: ToNucleotideValue, I: DoubleEndedIterator<Item = N>> DoubleEndedIterator
    for ReverseComplement<I>
{
    fn next_back(&mut self) -> Option<Self::Item> {
        self.0.next().map(|n| n.to_nucleotide().complement())
    }
}

impl<N: ToNucleotideValue, I: DoubleEndedIterator<Item = N> + ExactSizeIterator> ExactSizeIterator
    for ReverseComplement<I>
{
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
        let rc_frame1: Vec<_> = array.reverse_complement().frames()[1]
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

        let unchanged: Vec<_> = array
            .reverse_complement()
            .reverse_complement()
            .collect();
        assert_eq!(unchanged, array);
    }
}
