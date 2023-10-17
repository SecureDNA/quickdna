// Copyright 2021-2023 SecureDNA Stiftung (SecureDNA Foundation) <licensing@securedna.org>
// SPDX-License-Identifier: MIT OR Apache-2.0

use std::sync::Arc;

use smallvec::SmallVec;

use crate::{
    BaseSequence, DnaSequenceAmbiguous, DnaSequenceStrict, Nucleotide, NucleotideAmbiguous,
    NucleotideLike,
};

/// Iterator of all unambiguous expansions of ambiguous DNA.
///
/// Expansions are returned in lexicographic order based on the ordering of [`Nucleotide`]
/// (not currently alphabetical).
#[derive(Clone)]
pub struct Expansions {
    ambiguities: SmallVec<[Ambiguity; 8]>,
    buf: Arc<[Nucleotide]>,
    began: bool,
}

#[derive(Clone)]
struct Ambiguity {
    index: usize,
    digit: u8,
    nucleotide: NucleotideAmbiguous,
}

/// Unambiguous DNA expansion produced by [`Expansions`] iterator.
#[derive(Clone, PartialEq, Eq)]
pub struct Expansion(Arc<[Nucleotide]>);

impl Expansions {
    // Construct new [`Expansions`] iterator
    pub fn new(dna: &[NucleotideAmbiguous]) -> Self {
        let ambiguities = dna
            .iter()
            .enumerate()
            .filter(|(_, nuc)| nuc.bits().count_ones() > 1)
            .map(|(index, &nucleotide)| Ambiguity {
                index,
                digit: 0,
                nucleotide,
            })
            .collect();
        let buf: SmallVec<[_; 64]> = dna
            .iter()
            .map(|nuc| *nuc.possibilities().first().unwrap())
            .collect();
        let buf = Arc::from(buf.as_slice());
        Expansions {
            ambiguities,
            buf,
            began: false,
        }
    }
}

impl Iterator for Expansions {
    type Item = Expansion;

    fn next(&mut self) -> Option<Self::Item> {
        if !self.began {
            self.began = true;
            return Some(Expansion(self.buf.clone()));
        }

        // Can't use Arc::make_mut because [T] is unsized
        let buf = match Arc::get_mut(&mut self.buf) {
            Some(buf) => buf,
            None => {
                self.buf = Arc::from(&*self.buf);
                Arc::get_mut(&mut self.buf).unwrap()
            }
        };

        for amb in self.ambiguities.iter_mut().rev() {
            amb.digit = (amb.digit + 1) % (amb.nucleotide.bits().count_ones() as u8);
            buf[amb.index] = amb.nucleotide.possibilities()[amb.digit as usize];
            if amb.digit > 0 {
                return Some(Expansion(self.buf.clone()));
            }
        }
        None
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let size = (|| {
            let mut size: usize = 0;
            for amb in &self.ambiguities {
                let num_digit_states = amb.nucleotide.bits().count_ones() as usize;
                let remaining = num_digit_states - (amb.digit as usize) - 1;
                size = size.checked_mul(num_digit_states)?.checked_add(remaining)?;
            }
            size.checked_add(!self.began as usize)
        })();
        (size.unwrap_or(usize::MAX), size)
    }
}

impl std::fmt::Debug for Expansions {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let dna = self.buf.iter().map(|&nuc| nuc.into()).collect();
        let mut dna = DnaSequenceAmbiguous::new(dna);
        for amb in &self.ambiguities {
            dna[amb.index] = amb.nucleotide;
        }
        f.debug_tuple("Expansions").field(&dna.to_string()).finish()
    }
}

impl Expansion {
    /// Produces a [`DnaSequenceStrict`] from [`Expansion`].
    ///
    /// This takes *O*(*N*) time.
    pub fn to_dna(&self) -> DnaSequenceStrict {
        DnaSequenceStrict::new(self.to_vec())
    }
}

impl From<Expansion> for DnaSequenceStrict {
    fn from(expansion: Expansion) -> Self {
        expansion.to_dna()
    }
}

impl std::ops::Deref for Expansion {
    type Target = [Nucleotide];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl AsRef<[Nucleotide]> for Expansion {
    fn as_ref(&self) -> &[Nucleotide] {
        self
    }
}

impl std::cmp::PartialEq<&[Nucleotide]> for Expansion {
    fn eq(&self, other: &&[Nucleotide]) -> bool {
        self.as_ref() == *other
    }
}

impl std::cmp::PartialEq<DnaSequenceStrict> for Expansion {
    fn eq(&self, other: &DnaSequenceStrict) -> bool {
        self.as_ref() == other.as_slice()
    }
}

impl std::fmt::Debug for Expansion {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        f.debug_tuple("Expansion")
            .field(&self.to_dna().to_string())
            .finish()
    }
}

#[cfg(test)]
mod test {
    use super::*;

    use crate::{BaseSequence, DnaSequenceAmbiguous};

    fn amb_dna(dna: &str) -> DnaSequenceAmbiguous {
        dna.parse().unwrap()
    }

    fn dna(dna: &str) -> DnaSequenceStrict {
        dna.parse().unwrap()
    }

    #[test]
    fn basic_end_to_end() {
        let src_dna = amb_dna("ATBCGYAC"); // AT{TCG}CG{TC}AC
        let mut expansions = src_dna.expansions();
        assert_eq!(expansions.size_hint(), (6, Some(6)));
        assert_eq!(expansions.next().unwrap(), dna("ATTCGTAC"));
        assert_eq!(expansions.size_hint(), (5, Some(5)));
        assert_eq!(expansions.next().unwrap(), dna("ATTCGCAC"));
        assert_eq!(expansions.size_hint(), (4, Some(4)));
        assert_eq!(expansions.next().unwrap(), dna("ATCCGTAC"));
        assert_eq!(expansions.size_hint(), (3, Some(3)));
        assert_eq!(expansions.next().unwrap(), dna("ATCCGCAC"));
        assert_eq!(expansions.size_hint(), (2, Some(2)));
        assert_eq!(expansions.next().unwrap(), dna("ATGCGTAC"));
        assert_eq!(expansions.size_hint(), (1, Some(1)));
        assert_eq!(expansions.next().unwrap(), dna("ATGCGCAC"));
        assert_eq!(expansions.size_hint(), (0, Some(0)));
        assert!(expansions.next().is_none());
    }

    #[test]
    fn verify_basic_counting() {
        let actual: Vec<_> = amb_dna("WWW").expansions().collect();
        let expected = ["AAA", "AAT", "ATA", "ATT", "TAA", "TAT", "TTA", "TTT"].map(dna);
        assert_eq!(actual, expected);
    }

    #[test]
    fn verify_size_hints() {
        let testcases = [
            ("", 1),
            ("A", 1),
            ("ATCG", 1),
            ("W", 2),
            ("B", 3),
            ("N", 4),
            ("MR", 4),
            ("YV", 6),
            ("DS", 6),
            ("KN", 8),
            ("NW", 8),
            ("MRY", 8),
            ("HB", 9),
            ("ACGTYTAGCNCGATVGCTA", 24),
        ];
        for (src_dna, mut expected_len) in testcases {
            let mut iter = amb_dna(src_dna).expansions();
            assert_eq!(
                iter.size_hint(),
                (expected_len, Some(expected_len)),
                "Wrong initial size_hint for {src_dna:?}"
            );
            while iter.next().is_some() {
                assert!(expected_len > 0, "Iter too long for {src_dna:?}");
                expected_len -= 1;
                assert_eq!(
                    iter.size_hint(),
                    (expected_len, Some(expected_len)),
                    "Wrong size_hint during iteration of {src_dna:?}"
                );
            }
            assert_eq!(expected_len, 0, "Iterator too short for {src_dna:?}");
        }
    }

    #[test]
    fn unambiguous_sequences_have_single_element() {
        let src_dna = amb_dna("ATCGATATCGCGAATTCCGG");
        let mut expansions = src_dna.expansions();
        assert_eq!(expansions.size_hint(), (1, Some(1)));
        assert_eq!(expansions.next().unwrap(), dna("ATCGATATCGCGAATTCCGG"));
        assert_eq!(expansions.size_hint(), (0, Some(0)));
        assert!(expansions.next().is_none());
    }

    #[test]
    fn size_hint_handles_overflow() {
        let src_dna = amb_dna("NNNNNNNNNNNNNNNNNNBBBBBBBBBBBBBBBBBY");
        assert_eq!(
            src_dna.expansions().size_hint(),
            (17748888853923495936, Some(17748888853923495936))
        );
        let src_dna = amb_dna("NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN");
        assert_eq!(src_dna.expansions().size_hint(), (usize::MAX, None));
    }

    #[test]
    fn works_even_if_expansions_are_held() {
        let src_dna = amb_dna("ATCGNGCTA");
        let expansions: Vec<_> = src_dna.expansions().collect();
        let expected = [
            dna("ATCGAGCTA"),
            dna("ATCGTGCTA"),
            dna("ATCGCGCTA"),
            dna("ATCGGGCTA"),
        ];
        assert_eq!(expansions, expected);
    }

    #[test]
    fn debug_expansion() {
        let expansion = Expansion(Arc::from(dna("ATCG").as_slice()));
        assert_eq!(format!("{expansion:?}"), "Expansion(\"ATCG\")");
    }

    #[test]
    fn debug_expansions() {
        let src_dna = amb_dna("ATNCG");
        let expansions = Expansions::new(src_dna.as_slice());
        assert_eq!(format!("{expansions:?}"), "Expansions(\"ATNCG\")");
    }
}
