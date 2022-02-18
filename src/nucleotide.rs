use std::fmt::{self, Write};

use crate::errors::TranslationError;

/// A DNA nucleotide, or the IUPAC ambiguity code 'N'.
///
/// Sorts in ATCGN order, not alphabetical.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, std::hash::Hash)]
#[repr(u8)]
pub enum Nucleotide {
    A = 0,
    T = 1,
    C = 2,
    G = 3,
    /// IUPAC ambiguity code
    N = 4,
}

impl Nucleotide {
    pub const NUCLEOTIDES: [Self; 5] = [Self::A, Self::T, Self::C, Self::G, Self::N];
    const N_NUCLEOTIDES: u8 = Self::NUCLEOTIDES.len() as u8;

    const fn ascii_pack_table() -> [u8; 128] {
        let mut pack_table = [255u8; 128];
        pack_table[b'a' as usize] = Self::A as u8;
        pack_table[b'A' as usize] = Self::A as u8;
        pack_table[b't' as usize] = Self::T as u8;
        pack_table[b'T' as usize] = Self::T as u8;
        pack_table[b'c' as usize] = Self::C as u8;
        pack_table[b'C' as usize] = Self::C as u8;
        pack_table[b'g' as usize] = Self::G as u8;
        pack_table[b'G' as usize] = Self::G as u8;
        pack_table[b'n' as usize] = Self::N as u8;
        pack_table[b'N' as usize] = Self::N as u8;
        pack_table
    }

    const ASCII_PACK_TABLE: [u8; 128] = Self::ascii_pack_table();

    const COMPLEMENT_TABLE: [Self; Self::N_NUCLEOTIDES as usize] =
        [Self::T, Self::A, Self::G, Self::C, Self::N];
    pub const fn complement(self) -> Self {
        Self::COMPLEMENT_TABLE[self as u8 as usize]
    }

    const ASCII_MAP: [u8; Self::N_NUCLEOTIDES as usize] = [b'A', b'T', b'C', b'G', b'N'];
    pub const fn to_ascii(self) -> u8 {
        Self::ASCII_MAP[self as usize]
    }
}

impl TryFrom<u8> for Nucleotide {
    type Error = TranslationError;

    #[inline(always)]
    fn try_from(u: u8) -> Result<Self, Self::Error> {
        if std::intrinsics::likely(u < 128) {
            let v = Self::ASCII_PACK_TABLE[u as usize];
            if std::intrinsics::likely(v < Self::N_NUCLEOTIDES) {
                // SAFETY: there are only X variants, with assigned numbers, so 0..=(X - 1)
                //         are valid reprs of this type
                Ok(unsafe { std::mem::transmute(v) })
            } else {
                Err(TranslationError::BadNucleotide(u.into()))
            }
        } else {
            Err(TranslationError::NonAsciiByte(u))
        }
    }
}

impl From<Nucleotide> for u8 {
    fn from(n: Nucleotide) -> Self {
        n.to_ascii()
    }
}

impl From<Nucleotide> for char {
    fn from(n: Nucleotide) -> Self {
        n.to_ascii() as char
    }
}

impl fmt::Display for Nucleotide {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_char((*self).into())
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, std::hash::Hash)]
pub struct Codon(pub [Nucleotide; 3]);

impl Codon {
    /// How many of the specified nucleotide this codon includes
    pub fn count(&self, n: Nucleotide) -> usize {
        let mut count = 0;
        if self.0[0] == n {
            count += 1
        }
        if self.0[1] == n {
            count += 1
        }
        if self.0[2] == n {
            count += 1
        }
        count
    }

    /// Whether the codon includes `N`
    pub fn is_ambiguous(&self) -> bool {
        self.count(Nucleotide::N) > 0
    }

    /// Returns an iterator of all codons, including ambiguous codons with `N`
    pub fn all_codons() -> impl Iterator<Item = Self> {
        const N: [Nucleotide; 5] = Nucleotide::NUCLEOTIDES;
        N.iter().flat_map(move |&a| {
            N.iter()
                .flat_map(move |&b| N.iter().map(move |&c| Codon([a, b, c])))
        })
    }
}

impl TryFrom<[u8; 3]> for Codon {
    type Error = TranslationError;

    fn try_from(value: [u8; 3]) -> Result<Self, Self::Error> {
        Ok(Self([
            Nucleotide::try_from(value[0])?,
            Nucleotide::try_from(value[1])?,
            Nucleotide::try_from(value[2])?,
        ]))
    }
}

impl From<Codon> for [Nucleotide; 3] {
    fn from(c: Codon) -> Self {
        c.0
    }
}

impl fmt::Display for Codon {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}{}{}", self.0[0], self.0[1], self.0[2])
    }
}
