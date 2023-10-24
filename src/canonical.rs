// Copyright 2021-2023 SecureDNA Stiftung (SecureDNA Foundation) <licensing@securedna.org>
// SPDX-License-Identifier: MIT OR Apache-2.0

use std::cmp::Ordering;

use crate::Nucleotide;

/// Permute bases (and maybe reverse sequence) to produce lexical-minimum substitution of DNA.
///
/// This is similar to [`ForwardCanonical`] except in addition to remapping bases, it may also
/// reverse the original [`Nucleotide`] sequence, if doing so yields a lexically "earlier"
/// [`Nucleotide`] sequence. This means that `AATA` and `TCTT` both have the same canonical
/// sequence.
///
/// Thus, two [`Nucleotide`] sequence have the same canonical form if-and-only-if one is
/// isomorphic to the other (or its reverse). Canonicalization is idempotent.
#[derive(Clone, Debug)]
pub struct Canonical<I>(LexicalMin<ForwardCanonical<I>, ForwardCanonical<std::iter::Rev<I>>>);

impl<I> Canonical<I>
where
    I: Iterator<Item = Nucleotide>,
{
    /// Create iter of canonical substition for `iterable`
    pub fn new(iterable: impl IntoIterator<IntoIter = I>) -> Self
    where
        I: DoubleEndedIterator<Item = Nucleotide> + Clone,
    {
        let iter = iterable.into_iter();
        let fw_canon = ForwardCanonical::new(iter.clone());
        let rev_canon = ForwardCanonical::new(iter.rev());
        Canonical(LexicalMin::new(fw_canon, rev_canon))
    }
}

impl<I> Iterator for Canonical<I>
where
    I: DoubleEndedIterator<Item = Nucleotide>,
{
    type Item = Nucleotide;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        self.0.next()
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.0.size_hint()
    }
}

impl<I> ExactSizeIterator for Canonical<I> where
    I: ExactSizeIterator + DoubleEndedIterator<Item = Nucleotide>
{
}

/// Permute bases to produce lexical-minimum substitution of DNA.
///
/// This returns a sequence of [`Nucleotide`]s that is:
/// * Isomorphic to the original sequence; that is, the bases can be remapped to convert between
///   the original and forward-canonical sequences. So `ATA` and `GAG` have the same
///   forward-canonical sequence but `TAA` and `GAG` do not.
/// * Lexically "before" all other [`Nucleotide`] sequences that are isomorphic to it. Note that
///   (at the time of this writing) the order of [`Nucleotide`]s is
///   [`A`](Nucleotide::A) [`T`](Nucleotide::T) [`C`](Nucleotide::C) [`G`](Nucleotide::G)
///   instead of alphabetical, so the forward-canonical sequence for `CATTAG` is `ATCCTG`, _not_
///   `ACGGCT`.
///
/// Thus, two sequences of [`Nucleotide`]s produce the same forward-canonical sequence
/// if-and-only-if they are isomorphic. Forward-canonicalization is idempotent.
#[derive(Clone, Debug)]
pub struct ForwardCanonical<I> {
    inner: I,
    permutation: [Option<Nucleotide>; 4],
    unmapped: std::slice::Iter<'static, Nucleotide>,
}

impl<I> ForwardCanonical<I>
where
    I: Iterator<Item = Nucleotide>,
{
    /// Create iter of forward-canonical substition for `iterable`
    pub fn new(iterable: impl IntoIterator<IntoIter = I>) -> Self
    where
        I: Iterator<Item = Nucleotide>, // Might as well catch type errors early.
    {
        Self {
            inner: iterable.into_iter(),
            permutation: [None; 4],
            unmapped: Nucleotide::ALL.iter(),
        }
    }
}

impl<I> Iterator for ForwardCanonical<I>
where
    I: Iterator<Item = Nucleotide>,
{
    type Item = Nucleotide;

    fn next(&mut self) -> Option<Self::Item> {
        let idx = match self.inner.next()? {
            Nucleotide::A => 0,
            Nucleotide::T => 1,
            Nucleotide::C => 2,
            Nucleotide::G => 3,
        };
        let nuc = *self.permutation[idx]
            .get_or_insert_with(|| self.unmapped.next().copied().unwrap_or(Nucleotide::A));
        Some(nuc)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.inner.size_hint()
    }
}

impl<I: ExactSizeIterator<Item = Nucleotide>> ExactSizeIterator for ForwardCanonical<I> {}

// Given two sequences, returns whichever one is lexically less than the other.
// This is like an allocation-free equivalent of:
//     Vec::from_iter(iter1).min(iter2.collect()).into_iter()
// For example:
//     let x = [2, 1, 3];
//     let y = [2, 2, 2];
//     let lmin = LexicalMin::new(x.iter(), y.iter());
//     assert!(x < y);
//     assert!(lmin.eq(x.iter())); // because x < y
#[derive(Clone, Debug)]
struct LexicalMin<I1, I2> {
    iter1: I1,
    order: Ordering,
    iter2: I2,
}

impl<I1, I2> LexicalMin<I1, I2> {
    fn new(iter1: I1, iter2: I2) -> Self {
        Self {
            iter1,
            order: Ordering::Equal,
            iter2,
        }
    }
}

impl<I1, I2> Iterator for LexicalMin<I1, I2>
where
    I1: Iterator,
    I2: Iterator<Item = I1::Item>,
    I1::Item: Ord,
{
    type Item = I1::Item;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        match self.order {
            Ordering::Less => self.iter1.next(),
            Ordering::Equal => {
                let item1 = self.iter1.next();
                let item2 = self.iter2.next();
                self.order = item1.cmp(&item2);
                if self.order.is_lt() {
                    item1
                } else {
                    item2
                }
            }
            Ordering::Greater => self.iter2.next(),
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let (min1, max1) = self.iter1.size_hint();
        let (min2, max2) = self.iter2.size_hint();
        let min = min1.min(min2);
        let max = max1.zip(max2).map(|(max1, max2)| max1.max(max2));
        (min, max)
    }
}

#[cfg(test)]
mod test {
    use super::*;

    use quickcheck::quickcheck;

    use crate::{BaseSequence, DnaSequenceStrict};

    fn fw_canon(src_dna: &str) -> String {
        let dna: DnaSequenceStrict = src_dna.parse().unwrap();
        let canonical = ForwardCanonical::new(dna.iter()).collect();
        DnaSequenceStrict::new(canonical).to_string()
    }

    fn canon(src_dna: &str) -> String {
        let dna: DnaSequenceStrict = src_dna.parse().unwrap();
        dna.canonical().to_string()
    }

    #[test]
    fn sanity_check_forward_canonicalization() {
        assert_eq!(
            fw_canon("TGCGAGTGTAGCGAGATGTAGCGTAGAGTCTGAGATGCAGTA"),
            "ATCTGTATAGTCTGTGATAGTCTAGTGTACATGTGATCGTAG"
        );
    }

    #[test]
    fn sanity_check_canonicalization() {
        assert_eq!(
            canon("TGCGAGTGTAGCGAGATGTAGCGTAGAGTCTGAGATGCAGTA"),
            "ATCAGCTACACTGTCACATCGCATCTACACGCATCTCACGCT"
        );
    }

    #[test]
    fn exhaustively_check_forward_canonicalization_of_small_dna() {
        assert_eq!(fw_canon(""), "");

        assert_eq!(fw_canon("A"), "A");
        assert_eq!(fw_canon("T"), "A");
        assert_eq!(fw_canon("C"), "A");
        assert_eq!(fw_canon("G"), "A");

        assert_eq!(fw_canon("AA"), "AA");
        assert_eq!(fw_canon("AT"), "AT");
        assert_eq!(fw_canon("AC"), "AT");
        assert_eq!(fw_canon("AG"), "AT");
        assert_eq!(fw_canon("TA"), "AT");
        assert_eq!(fw_canon("TT"), "AA");
        assert_eq!(fw_canon("TC"), "AT");
        assert_eq!(fw_canon("TG"), "AT");
        assert_eq!(fw_canon("CA"), "AT");
        assert_eq!(fw_canon("CT"), "AT");
        assert_eq!(fw_canon("CC"), "AA");
        assert_eq!(fw_canon("CG"), "AT");
        assert_eq!(fw_canon("GA"), "AT");
        assert_eq!(fw_canon("GT"), "AT");
        assert_eq!(fw_canon("GC"), "AT");
        assert_eq!(fw_canon("GG"), "AA");

        assert_eq!(fw_canon("AAA"), "AAA");
        assert_eq!(fw_canon("AAT"), "AAT");
        assert_eq!(fw_canon("AAC"), "AAT");
        assert_eq!(fw_canon("AAG"), "AAT");
        assert_eq!(fw_canon("ATA"), "ATA");
        assert_eq!(fw_canon("ATT"), "ATT");
        assert_eq!(fw_canon("ATC"), "ATC");
        assert_eq!(fw_canon("ATG"), "ATC");
        assert_eq!(fw_canon("ACA"), "ATA");
        assert_eq!(fw_canon("ACT"), "ATC");
        assert_eq!(fw_canon("ACC"), "ATT");
        assert_eq!(fw_canon("ACG"), "ATC");
        assert_eq!(fw_canon("AGA"), "ATA");
        assert_eq!(fw_canon("AGT"), "ATC");
        assert_eq!(fw_canon("AGC"), "ATC");
        assert_eq!(fw_canon("AGG"), "ATT");

        assert_eq!(fw_canon("TAA"), "ATT");
        assert_eq!(fw_canon("TAT"), "ATA");
        assert_eq!(fw_canon("TAC"), "ATC");
        assert_eq!(fw_canon("TAG"), "ATC");
        assert_eq!(fw_canon("TTA"), "AAT");
        assert_eq!(fw_canon("TTT"), "AAA");
        assert_eq!(fw_canon("TTC"), "AAT");
        assert_eq!(fw_canon("TTG"), "AAT");
        assert_eq!(fw_canon("TCA"), "ATC");
        assert_eq!(fw_canon("TCT"), "ATA");
        assert_eq!(fw_canon("TCC"), "ATT");
        assert_eq!(fw_canon("TCG"), "ATC");
        assert_eq!(fw_canon("TGA"), "ATC");
        assert_eq!(fw_canon("TGT"), "ATA");
        assert_eq!(fw_canon("TGC"), "ATC");
        assert_eq!(fw_canon("TGG"), "ATT");

        assert_eq!(fw_canon("CAA"), "ATT");
        assert_eq!(fw_canon("CAT"), "ATC");
        assert_eq!(fw_canon("CAC"), "ATA");
        assert_eq!(fw_canon("CAG"), "ATC");
        assert_eq!(fw_canon("CTA"), "ATC");
        assert_eq!(fw_canon("CTT"), "ATT");
        assert_eq!(fw_canon("CTC"), "ATA");
        assert_eq!(fw_canon("CTG"), "ATC");
        assert_eq!(fw_canon("CCA"), "AAT");
        assert_eq!(fw_canon("CCT"), "AAT");
        assert_eq!(fw_canon("CCC"), "AAA");
        assert_eq!(fw_canon("CCG"), "AAT");
        assert_eq!(fw_canon("CGA"), "ATC");
        assert_eq!(fw_canon("CGT"), "ATC");
        assert_eq!(fw_canon("CGC"), "ATA");
        assert_eq!(fw_canon("CGG"), "ATT");

        assert_eq!(fw_canon("GAA"), "ATT");
        assert_eq!(fw_canon("GAT"), "ATC");
        assert_eq!(fw_canon("GAC"), "ATC");
        assert_eq!(fw_canon("GAG"), "ATA");
        assert_eq!(fw_canon("GTA"), "ATC");
        assert_eq!(fw_canon("GTT"), "ATT");
        assert_eq!(fw_canon("GTC"), "ATC");
        assert_eq!(fw_canon("GTG"), "ATA");
        assert_eq!(fw_canon("GCA"), "ATC");
        assert_eq!(fw_canon("GCT"), "ATC");
        assert_eq!(fw_canon("GCC"), "ATT");
        assert_eq!(fw_canon("GCG"), "ATA");
        assert_eq!(fw_canon("GGA"), "AAT");
        assert_eq!(fw_canon("GGT"), "AAT");
        assert_eq!(fw_canon("GGC"), "AAT");
        assert_eq!(fw_canon("GGG"), "AAA");
    }

    #[test]
    fn exhaustively_check_canonicalization_of_small_dna() {
        assert_eq!(canon(""), "");

        assert_eq!(canon("A"), "A");
        assert_eq!(canon("T"), "A");
        assert_eq!(canon("C"), "A");
        assert_eq!(canon("G"), "A");

        assert_eq!(canon("AA"), "AA");
        assert_eq!(canon("AT"), "AT");
        assert_eq!(canon("AC"), "AT");
        assert_eq!(canon("AG"), "AT");
        assert_eq!(canon("TA"), "AT");
        assert_eq!(canon("TT"), "AA");
        assert_eq!(canon("TC"), "AT");
        assert_eq!(canon("TG"), "AT");
        assert_eq!(canon("CA"), "AT");
        assert_eq!(canon("CT"), "AT");
        assert_eq!(canon("CC"), "AA");
        assert_eq!(canon("CG"), "AT");
        assert_eq!(canon("GA"), "AT");
        assert_eq!(canon("GT"), "AT");
        assert_eq!(canon("GC"), "AT");
        assert_eq!(canon("GG"), "AA");

        assert_eq!(canon("AAA"), "AAA");
        assert_eq!(canon("AAT"), "AAT");
        assert_eq!(canon("AAC"), "AAT");
        assert_eq!(canon("AAG"), "AAT");
        assert_eq!(canon("ATA"), "ATA");
        assert_eq!(canon("ATT"), "AAT");
        assert_eq!(canon("ATC"), "ATC");
        assert_eq!(canon("ATG"), "ATC");
        assert_eq!(canon("ACA"), "ATA");
        assert_eq!(canon("ACT"), "ATC");
        assert_eq!(canon("ACC"), "AAT");
        assert_eq!(canon("ACG"), "ATC");
        assert_eq!(canon("AGA"), "ATA");
        assert_eq!(canon("AGT"), "ATC");
        assert_eq!(canon("AGC"), "ATC");
        assert_eq!(canon("AGG"), "AAT");

        assert_eq!(canon("TAA"), "AAT");
        assert_eq!(canon("TAT"), "ATA");
        assert_eq!(canon("TAC"), "ATC");
        assert_eq!(canon("TAG"), "ATC");
        assert_eq!(canon("TTA"), "AAT");
        assert_eq!(canon("TTT"), "AAA");
        assert_eq!(canon("TTC"), "AAT");
        assert_eq!(canon("TTG"), "AAT");
        assert_eq!(canon("TCA"), "ATC");
        assert_eq!(canon("TCT"), "ATA");
        assert_eq!(canon("TCC"), "AAT");
        assert_eq!(canon("TCG"), "ATC");
        assert_eq!(canon("TGA"), "ATC");
        assert_eq!(canon("TGT"), "ATA");
        assert_eq!(canon("TGC"), "ATC");
        assert_eq!(canon("TGG"), "AAT");

        assert_eq!(canon("CAA"), "AAT");
        assert_eq!(canon("CAT"), "ATC");
        assert_eq!(canon("CAC"), "ATA");
        assert_eq!(canon("CAG"), "ATC");
        assert_eq!(canon("CTA"), "ATC");
        assert_eq!(canon("CTT"), "AAT");
        assert_eq!(canon("CTC"), "ATA");
        assert_eq!(canon("CTG"), "ATC");
        assert_eq!(canon("CCA"), "AAT");
        assert_eq!(canon("CCT"), "AAT");
        assert_eq!(canon("CCC"), "AAA");
        assert_eq!(canon("CCG"), "AAT");
        assert_eq!(canon("CGA"), "ATC");
        assert_eq!(canon("CGT"), "ATC");
        assert_eq!(canon("CGC"), "ATA");
        assert_eq!(canon("CGG"), "AAT");

        assert_eq!(canon("GAA"), "AAT");
        assert_eq!(canon("GAT"), "ATC");
        assert_eq!(canon("GAC"), "ATC");
        assert_eq!(canon("GAG"), "ATA");
        assert_eq!(canon("GTA"), "ATC");
        assert_eq!(canon("GTT"), "AAT");
        assert_eq!(canon("GTC"), "ATC");
        assert_eq!(canon("GTG"), "ATA");
        assert_eq!(canon("GCA"), "ATC");
        assert_eq!(canon("GCT"), "ATC");
        assert_eq!(canon("GCC"), "AAT");
        assert_eq!(canon("GCG"), "ATA");
        assert_eq!(canon("GGA"), "AAT");
        assert_eq!(canon("GGT"), "AAT");
        assert_eq!(canon("GGC"), "AAT");
        assert_eq!(canon("GGG"), "AAA");
    }

    #[test]
    fn canonicalization_selects_reverse_if_that_is_lexically_less() {
        assert_eq!(canon("TTGT"), "AATA");
        assert_eq!(canon("TGTT"), "AATA");

        assert_eq!(canon("ATCGCCAT"), "ATCCGCAT");
    }

    quickcheck! {
        fn forward_canonicalization_is_idempotent(dna: DnaSequenceStrict) -> bool {
            let canonical = ForwardCanonical::new(dna.as_slice().iter().copied());
            let canonical2 = ForwardCanonical::new(canonical.clone());
            canonical2.eq(canonical)
        }

        fn canonicalization_is_idempotent(dna: DnaSequenceStrict) -> bool {
            // Need full allocation because canonical needs a double-ended iter, but isn't one itself.
            let canonical = dna.canonical();
            let canonical2 = canonical.canonical();
            canonical2 == canonical
        }

        fn canonicalization_is_unaffected_by_reverse_complement(dna: DnaSequenceStrict) -> bool {
            dna.canonical() == dna.reverse_complement().canonical()
        }

        fn lexical_min_is_equivalent_to_vec_min(vec1: Vec<Nucleotide>, vec2: Vec<Nucleotide>) -> bool {
            let lmin = LexicalMin::new(vec1.iter(), vec2.iter());
            let vmin = vec1.clone().min(vec2.clone());
            lmin.eq(vmin.iter())
        }
    }
}
