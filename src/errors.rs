use thiserror::Error;

#[derive(Debug, Clone, Error)]
pub enum TranslationError {
    #[error("non-ascii byte: {:x?}", .0)]
    NonAsciiByte(u8),
    #[error("non-ascii char: {:x?}", .0)]
    NonAsciiChar(char),
    #[error("bad nucleotide: {:?}", .0)]
    BadNucleotide(char),
    #[error("bad amino acid: {:?}", .0)]
    BadAminoAcid(char),
    #[error("unexpected ambiguous nucleotide: {:?}", .0)]
    UnexpectedAmbiguousNucleotide(char),
    #[error("not a ncbi translation table: {}", .0)]
    BadTranslationTable(u8),
}
