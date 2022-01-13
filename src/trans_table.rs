use std::fmt;

use thiserror::Error;

#[derive(Debug, Clone, Error)]
pub enum TranslationError {
    #[error("non-ascii character: {:x?}", .0)]
    NonAscii(u8),
    #[error("bad nucleotide: {:?}", .0)]
    BadNucleotide(char),
    #[error("not a ncbi translation table: {}", .0)]
    BadTranslationTable(u8),
}

#[cold]
fn cold_error(error: TranslationError) -> TranslationError {
    error
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, std::hash::Hash)]
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
    pub fn complement(self) -> Self {
        Self::COMPLEMENT_TABLE[self as u8 as usize]
    }

    const ASCII_MAP: [u8; Self::N_NUCLEOTIDES as usize] = [b'A', b'T', b'C', b'G', b'N'];
    pub fn to_ascii(self) -> u8 {
        Self::ASCII_MAP[self as usize]
    }
}

impl TryFrom<u8> for Nucleotide {
    type Error = TranslationError;

    #[inline(always)]
    fn try_from(u: u8) -> Result<Self, Self::Error> {
        if u < 128 {
            let v = Self::ASCII_PACK_TABLE[u as usize];
            if v < Self::N_NUCLEOTIDES {
                // SAFETY: there are only X variants, with assigned numbers, so 0..=(X - 1)
                //         are valid reprs of this type
                Ok(unsafe { std::mem::transmute(v) })
            } else {
                Err(cold_error(TranslationError::BadNucleotide(u.into())))
            }
        } else {
            Err(cold_error(TranslationError::NonAscii(u)))
        }
    }
}

impl fmt::Display for Nucleotide {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::A => write!(f, "A"),
            Self::T => write!(f, "T"),
            Self::C => write!(f, "C"),
            Self::G => write!(f, "G"),
            Self::N => write!(f, "N"),
        }
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

#[repr(transparent)]
pub struct CodonIdx(usize);

impl From<[Nucleotide; 3]> for CodonIdx {
    fn from(value: [Nucleotide; 3]) -> Self {
        Self((value[0] as usize) << 6 | (value[1] as usize) << 3 | (value[2] as usize))
    }
}

impl From<Codon> for CodonIdx {
    fn from(c: Codon) -> Self {
        c.0.into()
    }
}

impl From<CodonIdx> for usize {
    fn from(c: CodonIdx) -> Self {
        c.0
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TranslationTable {
    Ncbi1,
    Ncbi2,
    Ncbi3,
    Ncbi4,
    Ncbi5,
    Ncbi6,
    Ncbi7,
    Ncbi8,
    Ncbi9,
    Ncbi10,
    Ncbi11,
    Ncbi12,
    Ncbi13,
    Ncbi14,
    Ncbi15,
    Ncbi16,
    // tables 17-20 are not assigned
    Ncbi21,
    Ncbi22,
    Ncbi23,
    Ncbi24,
    Ncbi25,
    Ncbi26,
    Ncbi27,
    Ncbi28,
    Ncbi29,
    Ncbi30,
    Ncbi31,
    Ncbi32,
    Ncbi33,
}

impl TranslationTable {
    /// There are really only 125 possible codons (len(ATCGN)^3), but since codons take up
    /// 3 bits, the maximum codon value is 0b100100100, so we need some holes in the table.
    pub const CODONS_PER_TABLE: usize = 293;
    // Number of NCBI translation tables (they go up to 33, but there's gaps in the numbering)
    pub const N_TRANS_TABLES: usize = 27;
    pub const LOOKUP_SIZE: usize = Self::CODONS_PER_TABLE * Self::N_TRANS_TABLES;
    /// Generated by bin/gen_table.rs, layout is all codons for table 1, then all codons for table 2, etc.
    const TRANSLATION_TABLES: &'static [u8; Self::LOOKUP_SIZE] = include_bytes!("tables.dat");

    fn table_index(self) -> usize {
        match self {
            // table 8 is an alias for table 1
            Self::Ncbi1 | Self::Ncbi8 => 0,
            Self::Ncbi2 => 1,
            Self::Ncbi3 => 2,
            // table 7 is identical to table 4
            Self::Ncbi4 | Self::Ncbi7 => 3,
            Self::Ncbi5 => 4,
            Self::Ncbi6 => 5,
            Self::Ncbi9 => 6,
            Self::Ncbi10 => 7,
            Self::Ncbi11 => 8,
            Self::Ncbi12 => 9,
            Self::Ncbi13 => 10,
            Self::Ncbi14 => 11,
            Self::Ncbi15 => 12,
            Self::Ncbi16 => 13,
            Self::Ncbi21 => 14,
            Self::Ncbi22 => 15,
            Self::Ncbi23 => 16,
            Self::Ncbi24 => 17,
            Self::Ncbi25 => 18,
            Self::Ncbi26 => 19,
            Self::Ncbi27 => 20,
            Self::Ncbi28 => 21,
            Self::Ncbi29 => 22,
            Self::Ncbi30 => 23,
            Self::Ncbi31 => 24,
            Self::Ncbi32 => 25,
            Self::Ncbi33 => 26,
        }
    }

    pub fn translate_dna(self, dna: &[u8]) -> Result<Vec<u8>, TranslationError> {
        if dna.is_empty() {
            return Ok(Vec::new());
        }

        let table_idx = self.table_index();

        let mut result = Vec::with_capacity(dna.len() / 3);

        // this will truncate any trailing non-multiple-of-3 chunk
        // biopython also truncates, but warns -- generally I don't think we care,
        // so I just made it silently truncate
        for &chunk in dna.array_chunks::<3>() {
            let a = chunk[0].try_into()?;
            let b = chunk[1].try_into()?;
            let c = chunk[2].try_into()?;
            let codon_idx = CodonIdx::from([a, b, c]);
            result.push(
                Self::TRANSLATION_TABLES
                    [table_idx * TranslationTable::CODONS_PER_TABLE + codon_idx.0],
            );
        }

        Ok(result)
    }
}

impl TryFrom<u8> for TranslationTable {
    type Error = TranslationError;

    fn try_from(value: u8) -> Result<Self, Self::Error> {
        match value {
            1 => Ok(Self::Ncbi1),
            2 => Ok(Self::Ncbi2),
            3 => Ok(Self::Ncbi3),
            4 => Ok(Self::Ncbi4),
            5 => Ok(Self::Ncbi5),
            6 => Ok(Self::Ncbi6),
            7 => Ok(Self::Ncbi7),
            8 => Ok(Self::Ncbi8),
            9 => Ok(Self::Ncbi9),
            10 => Ok(Self::Ncbi10),
            11 => Ok(Self::Ncbi11),
            12 => Ok(Self::Ncbi12),
            13 => Ok(Self::Ncbi13),
            14 => Ok(Self::Ncbi14),
            15 => Ok(Self::Ncbi15),
            16 => Ok(Self::Ncbi16),
            21 => Ok(Self::Ncbi21),
            22 => Ok(Self::Ncbi22),
            23 => Ok(Self::Ncbi23),
            24 => Ok(Self::Ncbi24),
            25 => Ok(Self::Ncbi25),
            26 => Ok(Self::Ncbi26),
            27 => Ok(Self::Ncbi27),
            28 => Ok(Self::Ncbi28),
            29 => Ok(Self::Ncbi29),
            30 => Ok(Self::Ncbi30),
            31 => Ok(Self::Ncbi31),
            32 => Ok(Self::Ncbi32),
            33 => Ok(Self::Ncbi33),
            _ => Err(TranslationError::BadTranslationTable(value)),
        }
    }
}

pub fn translate(table: u8, dna: &[u8]) -> Result<Vec<u8>, TranslationError> {
    let table = TranslationTable::try_from(table)?;
    table.translate_dna(dna)
}

pub fn reverse_complement(dna: &[u8]) -> Result<Vec<u8>, TranslationError> {
    let mut v = vec![0u8; dna.len()];
    for (i, &b) in dna.iter().enumerate() {
        let n = Nucleotide::try_from(b)?;
        v[dna.len() - 1 - i] = n.complement().to_ascii();
    }
    Ok(v)
}
