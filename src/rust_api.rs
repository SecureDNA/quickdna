use core::fmt;
use std::fmt::Write;
use std::str::FromStr;

use smallvec::SmallVec;

use crate::errors::TranslationError;
pub use crate::nucleotide::{Codon, Nucleotide};
pub use crate::trans_table::TranslationTable;

use crate::trans_table::reverse_complement;

pub trait BaseSequence: std::marker::Sized {
    type Item: Into<u8> + Copy;

    fn as_slice(&self) -> &[Self::Item];
    fn as_mut_slice(&mut self) -> &mut [Self::Item];

    fn len(&self) -> usize {
        self.as_slice().len()
    }

    fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

macro_rules! impls {
    ($type:ty) => {
        impl From<$type> for String {
            fn from(seq: $type) -> Self {
                let mut str = String::with_capacity(seq.len());
                for i in seq.as_slice() {
                    let c = u8::from(*i) as char;
                    str.push(c);
                }
                str
            }
        }

        impl std::ops::Index<usize> for $type {
            type Output = <Self as BaseSequence>::Item;

            fn index(&self, index: usize) -> &Self::Output {
                &self.as_slice()[index]
            }
        }

        impl std::ops::IndexMut<usize> for $type {
            fn index_mut(&mut self, index: usize) -> &mut Self::Output {
                &mut self.as_mut_slice()[index]
            }
        }

        impl $type {
            pub fn iter(&self) -> impl Iterator<Item = <Self as BaseSequence>::Item> + '_ {
                self.as_slice().iter().copied()
            }
        }
    };
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, std::hash::Hash)]
pub struct ProteinSequence {
    amino_acids: Vec<u8>,
}

impl ProteinSequence {
    pub fn new(amino_acids: Vec<u8>) -> Self {
        Self { amino_acids }
    }
}

impl BaseSequence for ProteinSequence {
    type Item = u8;

    fn as_slice(&self) -> &[u8] {
        &self.amino_acids
    }

    fn as_mut_slice(&mut self) -> &mut [u8] {
        &mut self.amino_acids
    }
}

impls!(ProteinSequence);

impl fmt::Display for ProteinSequence {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let str = String::from_utf8_lossy(&self.amino_acids);
        f.write_str(&str)
    }
}

impl FromStr for ProteinSequence {
    type Err = TranslationError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if std::intrinsics::likely(s.is_ascii()) {
            Ok(Self {
                amino_acids: s.as_bytes().to_vec(),
            })
        } else {
            let first_non_ascii = s.chars().find(|c| !c.is_ascii()).unwrap();
            Err(TranslationError::NonAsciiChar(first_non_ascii))
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, std::hash::Hash)]
pub struct DnaSequence {
    dna: Vec<Nucleotide>,
}

impl DnaSequence {
    /// Construct a new DnaSequence from a Vec of nucleotides
    pub fn new(dna: Vec<Nucleotide>) -> Self {
        Self { dna }
    }

    /// Translate this DNA sequence into a protein sequence, using the specified
    /// translation table.
    pub fn translate(&self, table: TranslationTable) -> ProteinSequence {
        let amino_acids = table.translate_dna(&self.dna);
        ProteinSequence::new(amino_acids)
    }

    /// Translate this DNA sequence into up to 3 protein sequences, one for each possible
    /// reading frame on this sense.
    ///
    /// May return less than 3 proteins for too-short sequences.
    /// For example, a sequence of length 4 only has 2 reading frames,
    /// and a sequence of length 2 has none.
    pub fn translate_self_frames(&self, table: TranslationTable) -> SmallVec<[ProteinSequence; 3]> {
        let mut result = SmallVec::with_capacity(3);

        // avoid empty translations & multiple branches
        if std::intrinsics::likely(self.len() >= 5) {
            result.push(ProteinSequence {
                amino_acids: table.translate_dna(&self.dna[0..]),
            });
            result.push(ProteinSequence {
                amino_acids: table.translate_dna(&self.dna[1..]),
            });
            result.push(ProteinSequence {
                amino_acids: table.translate_dna(&self.dna[2..]),
            });
        } else if self.len() == 4 {
            result.push(ProteinSequence {
                amino_acids: table.translate_dna(&self.dna[0..]),
            });
            result.push(ProteinSequence {
                amino_acids: table.translate_dna(&self.dna[1..]),
            });
        } else if self.len() == 3 {
            result.push(ProteinSequence {
                amino_acids: table.translate_dna(&self.dna[0..]),
            });
        }

        result
    }

    /// Translate this DNA sequence into at most 6 protein sequences, one for each possible
    /// reading frame on this sense and the reverse complement.
    ///
    /// May return less than 6 proteins for too-short sequences.
    /// For example, a sequence of length 4 only has 2 reading frames,
    /// and a sequence of length 2 has none.
    pub fn translate_all_frames(&self, table: TranslationTable) -> SmallVec<[ProteinSequence; 6]> {
        let mut result = SmallVec::with_capacity(6);

        result.append(&mut self.translate_self_frames(table));
        result.append(&mut self.reverse_complement().translate_self_frames(table));

        result
    }

    /// Takes the reverse complement of a DNA sequence.
    pub fn reverse_complement(&self) -> Self {
        Self::new(reverse_complement(&self.dna))
    }
}

impl BaseSequence for DnaSequence {
    type Item = Nucleotide;

    fn as_slice(&self) -> &[Self::Item] {
        &self.dna
    }

    fn as_mut_slice(&mut self) -> &mut [Self::Item] {
        &mut self.dna
    }
}

impls!(DnaSequence);

impl fmt::Display for DnaSequence {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for &n in &self.dna {
            f.write_char(n.into())?;
        }
        Ok(())
    }
}

impl FromStr for DnaSequence {
    type Err = TranslationError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut vec = vec![Nucleotide::N; s.len()];
        for (idx, &b) in s.as_bytes().iter().enumerate() {
            let n = Nucleotide::try_from(b)?;
            vec[idx] = n;
        }
        Ok(Self::new(vec))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use smallvec::smallvec;
    use std::hash::Hasher;

    macro_rules! assert_eq_smallvec {
        ($size:expr, $left:expr, $right:expr) => {
            let sv: SmallVec<[_; $size]> = $right;
            assert_eq!($left, sv);
        };
    }

    fn dna(dna: &str) -> DnaSequence {
        DnaSequence::from_str(dna).unwrap()
    }

    fn protein(aa: &str) -> ProteinSequence {
        ProteinSequence::from_str(aa).unwrap()
    }

    fn hash(v: impl std::hash::Hash) -> u64 {
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        v.hash(&mut hasher);
        hasher.finish()
    }

    #[test]
    fn test_translate() {
        assert_eq!(
            dna("AAAGGGAAA").translate(TranslationTable::Ncbi1),
            protein("KGK")
        );
    }

    #[test]
    fn test_translate_self() {
        assert_eq_smallvec!(
            3,
            dna("AAAGGGAAA").translate_self_frames(TranslationTable::Ncbi1),
            smallvec![protein("KGK"), protein("KG"), protein("RE"),]
        );
    }

    #[test]
    fn test_short_translate_self() {
        assert_eq_smallvec!(
            3,
            dna("GGGG").translate_self_frames(TranslationTable::Ncbi1),
            smallvec![protein("G"), protein("G"),]
        );

        assert_eq_smallvec!(
            3,
            dna("GGG").translate_self_frames(TranslationTable::Ncbi1),
            smallvec![protein("G"),]
        );

        assert!(dna("GG")
            .translate_self_frames(TranslationTable::Ncbi1)
            .is_empty());

        assert!(dna("G")
            .translate_self_frames(TranslationTable::Ncbi1)
            .is_empty());

        assert!(dna("")
            .translate_self_frames(TranslationTable::Ncbi1)
            .is_empty());
    }

    #[test]
    fn test_translate_all() {
        assert_eq_smallvec!(
            6,
            dna("AAAGGGAAA").translate_all_frames(TranslationTable::Ncbi1),
            smallvec![
                protein("KGK"),
                protein("KG"),
                protein("RE"),
                protein("FPF"),
                protein("FP"),
                protein("SL"),
            ]
        );
    }

    #[test]
    fn test_short_translate_all() {
        assert_eq_smallvec!(
            6,
            dna("GGGG").translate_all_frames(TranslationTable::Ncbi1),
            smallvec![protein("G"), protein("G"), protein("P"), protein("P"),]
        );

        assert_eq_smallvec!(
            6,
            dna("GGG").translate_all_frames(TranslationTable::Ncbi1),
            smallvec![protein("G"), protein("P"),]
        );

        assert!(dna("GG")
            .translate_all_frames(TranslationTable::Ncbi1)
            .is_empty());

        assert!(dna("G")
            .translate_all_frames(TranslationTable::Ncbi1)
            .is_empty());

        assert!(dna("")
            .translate_all_frames(TranslationTable::Ncbi1)
            .is_empty());
    }

    #[test]
    fn test_equality() {
        let d1 = dna("aaa");
        let d2 = dna("aaa");
        let p1 = protein("aaa");
        let p2 = protein("aaa");

        assert_eq!(d1, d2);
        assert_eq!(p1, p2);
    }

    #[test]
    fn test_hash() {
        let d1 = dna("aaa");
        let d2 = dna("aaa");
        let p1 = protein("aaa");
        let p2 = protein("aaa");

        assert_eq!(hash(&d1), hash(&d2));
        assert!(hash(&d1) != hash(&p1));
        assert_eq!(hash(&p1), hash(&p2));
    }
}
