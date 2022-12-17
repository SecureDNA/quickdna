use crate::{
    errors::TranslationError,
    nucleotide::{Codon, CodonAmbiguous, NucleotideLike},
};

/// Identifies a translation table for turning codons into amino acids.
/// See: <https://en.wikipedia.org/wiki/List_of_genetic_codes>
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TranslationTable {
    /// The standard code
    Ncbi1,
    /// The vertebrate mitochondrial code
    Ncbi2,
    /// The yeast mitochondrial code
    Ncbi3,
    /// The mold, protozoan, and coelenterate mitochondrial code and the mycoplasma/spiroplasma code
    Ncbi4,
    /// The invertebrate mitochondrial code
    Ncbi5,
    /// The ciliate, dasycladacean and hexamita nuclear code
    Ncbi6,
    /// The kinetoplast code; cf. table 4.
    Ncbi7,
    /// Same as table 1.
    Ncbi8,
    /// The echinoderm and flatworm mitochondrial code
    Ncbi9,
    /// The euplotid nuclear code
    Ncbi10,
    /// The bacterial, archaeal and plant plastid code
    Ncbi11,
    /// The alternative yeast nuclear code
    Ncbi12,
    /// The ascidian mitochondrial code
    Ncbi13,
    /// The alternative flatworm mitochondrial code
    Ncbi14,
    /// The Blepharisma nuclear code
    Ncbi15,
    /// The chlorophycean mitochondrial code
    Ncbi16,

    // tables 17-20 are not assigned
    /// The trematode mitochondrial code
    Ncbi21,
    /// The Scenedesmus obliquus mitochondrial code
    Ncbi22,
    /// The Thraustochytrium mitochondrial code
    Ncbi23,
    /// The Pterobranchia mitochondrial code
    Ncbi24,
    /// The candidate division SR1 and gracilibacteria code
    Ncbi25,
    /// The Pachysolen tannophilus nuclear code
    Ncbi26,
    /// The karyorelict nuclear code
    Ncbi27,
    /// The Condylostoma nuclear code
    Ncbi28,
    /// The Mesodinium nuclear code
    Ncbi29,
    /// The Peritrich nuclear code
    Ncbi30,
    /// The Blastocrithidia nuclear code
    Ncbi31,
    /// The Balanophoraceae plastid code
    Ncbi32,
    /// The Cephalodiscidae mitochondrial code
    Ncbi33,
}

#[repr(transparent)]
pub struct CodonIdx(usize);

impl<T: NucleotideLike> From<[T; 3]> for CodonIdx {
    fn from(value: [T; 3]) -> Self {
        let v0 = value[0].bits() as usize;
        let v1 = value[1].bits() as usize;
        let v2 = value[2].bits() as usize;
        Self((v0 << 8) | (v1 << 4) | v2)
    }
}

impl From<Codon> for CodonIdx {
    fn from(c: Codon) -> Self {
        c.0.into()
    }
}

impl From<CodonAmbiguous> for CodonIdx {
    fn from(c: CodonAmbiguous) -> Self {
        c.0.into()
    }
}

impl From<CodonIdx> for usize {
    fn from(c: CodonIdx) -> Self {
        c.0
    }
}

impl TranslationTable {
    /// Each ambiguity code is represented by 4 bits, so there are (2^4)^3 codons per table.
    pub const CODONS_PER_TABLE: usize = 1 << 12;
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

    /// Convert this table to a callable that maps codons to amino acids
    ///
    /// Currently, amino acids are represented as [`u8`]s containing the ascii
    /// code for the corresponding letter abbreviation.
    ///
    /// # Examples
    ///
    /// ```
    /// use quickdna::{Nucleotide, NucleotideIter, TranslationTable};
    ///
    /// use Nucleotide::*;
    /// let dna = [A, T, C, G, A, T, C, G];
    ///
    /// let ncbi1 = TranslationTable::Ncbi1.to_fn();
    /// let aas = dna.iter().codons().map(ncbi1);
    /// assert!(aas.eq([b'I', b'D']));
    /// ```
    pub fn to_fn<N: NucleotideLike, C: Into<[N; 3]>>(self) -> impl Copy + Fn(C) -> u8 {
        let start = self.table_index() * Self::CODONS_PER_TABLE;
        let end = start + Self::CODONS_PER_TABLE;
        let table = &Self::TRANSLATION_TABLES[start..end];
        |codon| {
            let nucleotides: [N; 3] = codon.into();
            let CodonIdx(i) = nucleotides.into();
            table[i]
        }
    }

    pub fn translate_dna_bytes<T: NucleotideLike>(
        self,
        dna: &[u8],
    ) -> Result<Vec<u8>, TranslationError> {
        if dna.is_empty() {
            return Ok(Vec::new());
        }

        let table_idx = self.table_index();

        let mut result = Vec::with_capacity(dna.len() / 3);

        // this will truncate any trailing non-multiple-of-3 chunk
        // biopython also truncates, but warns -- generally I don't think we care,
        // so I just made it silently truncate
        for chunk in dna.chunks_exact(3) {
            let a: T = chunk[0].try_into()?;
            let b: T = chunk[1].try_into()?;
            let c: T = chunk[2].try_into()?;
            let codon_idx = CodonIdx::from([a, b, c]);
            result.push(
                Self::TRANSLATION_TABLES
                    [table_idx * TranslationTable::CODONS_PER_TABLE + usize::from(codon_idx)],
            );
        }

        Ok(result)
    }

    pub fn translate_dna<T: NucleotideLike>(self, dna: &[T]) -> Vec<u8> {
        if dna.is_empty() {
            return Vec::new();
        }

        let table_idx = self.table_index();

        let mut result = Vec::with_capacity(dna.len() / 3);

        // this will truncate any trailing non-multiple-of-3 chunk
        // biopython also truncates, but warns -- generally I don't think we care,
        // so I just made it silently truncate
        for chunk in dna.chunks_exact(3) {
            let sized_chunk: [T; 3] = [chunk[0], chunk[1], chunk[2]];
            let codon_idx = CodonIdx::from(sized_chunk);
            result.push(
                Self::TRANSLATION_TABLES
                    [table_idx * TranslationTable::CODONS_PER_TABLE + usize::from(codon_idx)],
            );
        }

        result
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

pub fn reverse_complement_bytes<T: NucleotideLike>(
    dna: &[u8],
) -> Result<Vec<u8>, TranslationError> {
    let mut v = vec![0u8; dna.len()];
    for (i, &b) in dna.iter().enumerate() {
        let n = T::try_from(b)?;
        v[dna.len() - 1 - i] = n.complement().to_ascii();
    }
    Ok(v)
}

// Perf: it looks like .collect() gets great codegen here, but not when dealing
// with Result as it would in `reverse_complement_bytes` above. Here, it beats
// writing something like `let mut v = vec![T::default(), dna.len()];` which
// wastes time doing a memset before filling the Vec.
pub fn reverse_complement<T: NucleotideLike>(dna: &[T]) -> Vec<T> {
    dna.iter().rev().map(|n| n.complement()).collect()
}
