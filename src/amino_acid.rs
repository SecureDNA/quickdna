use crate::errors::TranslationError;

// TODO: Should I assign concrete numbers to these to ensure they can safely be hashed directly?
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, std::hash::Hash)]
#[repr(u8)]
pub enum AminoAcid {
    A,
    C,
    D,
    E,
    F,
    G,
    H,
    I,
    K,
    L,
    M,
    N,
    P,
    Q,
    R,
    S,
    T,
    V,
    W,
    Y,
}

const fn ascii_to_amino_acid_table() -> [Option<AminoAcid>; 256] {
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

    insert!(b'a', AminoAcid::A); // Alanine
    insert!(b'c', AminoAcid::C); // Cysteine
    insert!(b'd', AminoAcid::D); // Aspartate
    insert!(b'e', AminoAcid::E); // Glutamate
    insert!(b'f', AminoAcid::F); // Phenylalanine
    insert!(b'g', AminoAcid::G); // Glycine
    insert!(b'h', AminoAcid::H); // Histidine
    insert!(b'i', AminoAcid::I); // Isoleucine
    insert!(b'k', AminoAcid::K); // Lysine
    insert!(b'l', AminoAcid::L); // Leucine
    insert!(b'm', AminoAcid::M); // Methionine
    insert!(b'n', AminoAcid::N); // Asparagine
    insert!(b'p', AminoAcid::P); // Proline
    insert!(b'q', AminoAcid::Q); // Glutamine
    insert!(b'r', AminoAcid::R); // Arginine
    insert!(b's', AminoAcid::S); // Serine
    insert!(b't', AminoAcid::T); // Threonine
    insert!(b'v', AminoAcid::V); // Valine
    insert!(b'w', AminoAcid::W); // Tryptophan
    insert!(b'y', AminoAcid::Y); // Tyrosine

    pack_table
}

// FIXME: Workaround for AminoAcid::try_from::<u8>() not being const
pub const ASCII_TO_AMINO_ACID: [Option<AminoAcid>; 256] = ascii_to_amino_acid_table();

impl AminoAcid {
    pub fn to_ascii(self) -> u8 {
        match self {
            Self::A => b'A',
            Self::C => b'C',
            Self::D => b'D',
            Self::E => b'E',
            Self::F => b'F',
            Self::G => b'G',
            Self::H => b'H',
            Self::I => b'I',
            Self::K => b'K',
            Self::L => b'L',
            Self::M => b'M',
            Self::N => b'N',
            Self::P => b'P',
            Self::Q => b'Q',
            Self::R => b'R',
            Self::S => b'S',
            Self::T => b'T',
            Self::V => b'V',
            Self::W => b'W',
            Self::Y => b'Y',
        }
    }
}

impl TryFrom<u8> for AminoAcid {
    type Error = TranslationError;

    #[inline(always)]
    fn try_from(u: u8) -> Result<Self, Self::Error> {
        if u >= 128 {
            return Err(TranslationError::NonAsciiByte(u));
        }

        match ASCII_TO_AMINO_ACID[u as usize] {
            Some(na) => Ok(na),
            None => Err(TranslationError::BadAminoAcid(u.into())),
        }
    }
}
