use std::fmt::{self, Write};

use crate::errors::TranslationError;

/// A DNA nucleotide.
///
/// Sorts in ATCG order, not alphabetical.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, std::hash::Hash)]
#[repr(u8)]
pub enum Nucleotide {
    A = 0b0001,
    T = 0b0010,
    C = 0b0100,
    G = 0b1000,
}

/// A DNA nucleotide, or an IUPAC ambiguity code representing a set of possible nucleotides.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, std::hash::Hash)]
#[repr(u8)]
pub enum NucleotideAmbiguous {
    A = Nucleotide::A as u8,
    T = Nucleotide::T as u8,
    C = Nucleotide::C as u8,
    G = Nucleotide::G as u8,

    W = Nucleotide::A as u8 | Nucleotide::T as u8,
    M = Nucleotide::A as u8 | Nucleotide::C as u8,
    R = Nucleotide::A as u8 | Nucleotide::G as u8, // purines
    Y = Nucleotide::T as u8 | Nucleotide::C as u8, // pyrimidines
    S = Nucleotide::C as u8 | Nucleotide::G as u8,
    K = Nucleotide::T as u8 | Nucleotide::G as u8,

    B = Nucleotide::T as u8 | Nucleotide::C as u8 | Nucleotide::G as u8, // not A
    V = Nucleotide::A as u8 | Nucleotide::C as u8 | Nucleotide::G as u8, // not T
    D = Nucleotide::A as u8 | Nucleotide::T as u8 | Nucleotide::G as u8, // not C
    H = Nucleotide::A as u8 | Nucleotide::T as u8 | Nucleotide::C as u8, // not G

    N = Nucleotide::A as u8 | Nucleotide::T as u8 | Nucleotide::C as u8 | Nucleotide::G as u8,
}

pub trait NucleotideLike:
    Copy + Eq + Into<u8> + Into<char> + TryFrom<u8, Error = TranslationError>
{
    fn complement(self) -> Self;
    fn bits(self) -> u8;
    fn to_ascii(self) -> u8;
    fn is_ambiguous(self) -> bool;
}

const fn ascii_to_nucleotide_table() -> [Option<NucleotideAmbiguous>; 256] {
    // PERF: This was previously a 128-byte array, with "high" bytes eliminated
    // by a separate check. But extending it to 256 elements makes Rust realize
    // it's safe to index with any u8-cast-to-usize and eliminate a bounds check
    // from the generated code.
    let mut pack_table = [None; 256];

    macro_rules! insert {
        ($chr:literal, $variant:expr) => {
            pack_table[$chr.to_ascii_uppercase() as usize] = Some($variant);
            pack_table[$chr.to_ascii_lowercase() as usize] = Some($variant);
        };
    }

    insert!(b'a', NucleotideAmbiguous::A);
    insert!(b't', NucleotideAmbiguous::T);
    insert!(b'c', NucleotideAmbiguous::C);
    insert!(b'g', NucleotideAmbiguous::G);

    // ambiguity codes
    insert!(b'n', NucleotideAmbiguous::N); // A, T, C, or G (complement N)
    insert!(b'm', NucleotideAmbiguous::M); // A or C (complement K)
    insert!(b'r', NucleotideAmbiguous::R); // A or G (complement Y)
    insert!(b'w', NucleotideAmbiguous::W); // A or T (complement W, itself)
    insert!(b's', NucleotideAmbiguous::S); // C or G (complement S, itself)
    insert!(b'y', NucleotideAmbiguous::Y); // C or T (complement R)
    insert!(b'k', NucleotideAmbiguous::K); // G or T (complement M)
    insert!(b'v', NucleotideAmbiguous::V); // A, C or G (complement B)
    insert!(b'h', NucleotideAmbiguous::H); // A, C or T (complement D)
    insert!(b'd', NucleotideAmbiguous::D); // A, G or T (complement H)
    insert!(b'b', NucleotideAmbiguous::B); // C, G or T (complement V)

    pack_table
}

const ASCII_TO_NUCLEOTIDE: [Option<NucleotideAmbiguous>; 256] = ascii_to_nucleotide_table();

impl Nucleotide {
    pub const ALL: [Self; 4] = [Self::A, Self::T, Self::C, Self::G];
    pub const PURINES: [Self; 2] = [Self::A, Self::G];
    pub const PYRIMIDINES: [Self; 2] = [Self::C, Self::T];
}

impl NucleotideLike for Nucleotide {
    fn complement(self) -> Self {
        match self {
            Self::A => Self::T,
            Self::T => Self::A,
            Self::C => Self::G,
            Self::G => Self::C,
        }
    }

    fn bits(self) -> u8 {
        self as u8
    }

    fn to_ascii(self) -> u8 {
        match self {
            Self::A => b'A',
            Self::T => b'T',
            Self::C => b'C',
            Self::G => b'G',
        }
    }

    fn is_ambiguous(self) -> bool {
        false
    }
}

impl NucleotideAmbiguous {
    pub const ALL: [Self; 15] = [
        Self::A,
        Self::T,
        Self::W,
        Self::C,
        Self::M,
        Self::Y,
        Self::H,
        Self::G,
        Self::R,
        Self::K,
        Self::D,
        Self::S,
        Self::V,
        Self::B,
        Self::N,
    ];

    pub const fn possibilities(self) -> &'static [Nucleotide] {
        match self {
            Self::A => &[Nucleotide::A],
            Self::T => &[Nucleotide::T],
            Self::C => &[Nucleotide::C],
            Self::G => &[Nucleotide::G],
            Self::W => &[Nucleotide::A, Nucleotide::T],
            Self::M => &[Nucleotide::A, Nucleotide::C],
            Self::R => &[Nucleotide::A, Nucleotide::G],
            Self::Y => &[Nucleotide::T, Nucleotide::C],
            Self::S => &[Nucleotide::C, Nucleotide::G],
            Self::K => &[Nucleotide::T, Nucleotide::G],
            Self::B => &[Nucleotide::T, Nucleotide::C, Nucleotide::G],
            Self::V => &[Nucleotide::A, Nucleotide::C, Nucleotide::G],
            Self::D => &[Nucleotide::A, Nucleotide::T, Nucleotide::G],
            Self::H => &[Nucleotide::A, Nucleotide::T, Nucleotide::C],
            Self::N => &[Nucleotide::A, Nucleotide::T, Nucleotide::C, Nucleotide::G],
        }
    }
}

impl NucleotideLike for NucleotideAmbiguous {
    fn complement(self) -> Self {
        match self {
            Self::A => Self::T,
            Self::T => Self::A,
            Self::W => Self::W,
            Self::C => Self::G,
            Self::M => Self::K,
            Self::Y => Self::R,
            Self::H => Self::D,
            Self::G => Self::C,
            Self::R => Self::Y,
            Self::K => Self::M,
            Self::D => Self::H,
            Self::S => Self::S,
            Self::V => Self::B,
            Self::B => Self::V,
            Self::N => Self::N,
        }
    }

    fn bits(self) -> u8 {
        self as u8
    }

    fn to_ascii(self) -> u8 {
        match self {
            Self::A => b'A',
            Self::T => b'T',
            Self::W => b'W',
            Self::C => b'C',
            Self::M => b'M',
            Self::Y => b'Y',
            Self::H => b'H',
            Self::G => b'G',
            Self::R => b'R',
            Self::K => b'K',
            Self::D => b'D',
            Self::S => b'S',
            Self::V => b'V',
            Self::B => b'B',
            Self::N => b'N',
        }
    }

    fn is_ambiguous(self) -> bool {
        (self as usize).count_ones() > 1
    }
}

impl From<Nucleotide> for NucleotideAmbiguous {
    #[inline(always)]
    fn from(value: Nucleotide) -> Self {
        match value {
            Nucleotide::A => Self::A,
            Nucleotide::T => Self::T,
            Nucleotide::C => Self::C,
            Nucleotide::G => Self::G,
        }
    }
}

impl TryFrom<NucleotideAmbiguous> for Nucleotide {
    type Error = TranslationError;

    #[inline(always)]
    fn try_from(value: NucleotideAmbiguous) -> Result<Self, Self::Error> {
        match value {
            NucleotideAmbiguous::A => Ok(Self::A),
            NucleotideAmbiguous::T => Ok(Self::T),
            NucleotideAmbiguous::C => Ok(Self::C),
            NucleotideAmbiguous::G => Ok(Self::G),
            other => Err(TranslationError::UnexpectedAmbiguousNucleotide(
                other.into(),
            )),
        }
    }
}

impl TryFrom<u8> for Nucleotide {
    type Error = TranslationError;

    #[inline(always)]
    fn try_from(u: u8) -> Result<Self, Self::Error> {
        if u >= 128 {
            return Err(TranslationError::NonAsciiByte(u));
        }

        match ASCII_TO_NUCLEOTIDE[u as usize] {
            Some(na) => Nucleotide::try_from(na),
            None => Err(TranslationError::BadNucleotide(u.into())),
        }
    }
}

impl TryFrom<u8> for NucleotideAmbiguous {
    type Error = TranslationError;

    #[inline(always)]
    fn try_from(u: u8) -> Result<Self, Self::Error> {
        if u >= 128 {
            return Err(TranslationError::NonAsciiByte(u));
        }

        match ASCII_TO_NUCLEOTIDE[u as usize] {
            Some(na) => Ok(na),
            None => Err(TranslationError::BadNucleotide(u.into())),
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

impl From<NucleotideAmbiguous> for u8 {
    fn from(n: NucleotideAmbiguous) -> Self {
        n.to_ascii()
    }
}

impl From<NucleotideAmbiguous> for char {
    fn from(n: NucleotideAmbiguous) -> Self {
        n.to_ascii() as char
    }
}

impl fmt::Display for Nucleotide {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_char((*self).into())
    }
}

impl fmt::Display for NucleotideAmbiguous {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_char((*self).into())
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, std::hash::Hash)]
pub struct Codon(pub [Nucleotide; 3]);

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

#[derive(Debug, Clone, Copy, PartialEq, Eq, std::hash::Hash)]
pub struct CodonAmbiguous(pub [NucleotideAmbiguous; 3]);

impl TryFrom<[u8; 3]> for CodonAmbiguous {
    type Error = TranslationError;

    fn try_from(value: [u8; 3]) -> Result<Self, Self::Error> {
        Ok(Self([
            NucleotideAmbiguous::try_from(value[0])?,
            NucleotideAmbiguous::try_from(value[1])?,
            NucleotideAmbiguous::try_from(value[2])?,
        ]))
    }
}

impl From<CodonAmbiguous> for [NucleotideAmbiguous; 3] {
    fn from(c: CodonAmbiguous) -> Self {
        c.0
    }
}

impl fmt::Display for CodonAmbiguous {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}{}{}", self.0[0], self.0[1], self.0[2])
    }
}

impl CodonAmbiguous {
    pub fn possibilities(&self) -> impl Iterator<Item = Codon> + '_ {
        self.0[0].possibilities().iter().flat_map(move |&a| {
            self.0[1].possibilities().iter().flat_map(move |&b| {
                self.0[2]
                    .possibilities()
                    .iter()
                    .map(move |&c| Codon([a, b, c]))
            })
        })
    }
}
