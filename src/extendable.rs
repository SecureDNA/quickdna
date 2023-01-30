/// A trait used by the FASTA parser: `T: Extendable` is the type of the
/// contents of the FASTA file (String or DnaSequence or ProteinSequence).
/// This trait lets us generically concatenate parsed content lines.
pub trait Extendable {
    fn empty() -> Self;
    fn is_empty(&self) -> bool;
    fn extend(&mut self, other: Self) -> ();
}

impl Extendable for String {
    fn empty() -> Self {
        "".to_string()
    }

    /// When parsing a FASTA file, we consider whitespace-only lines to be "empty".
    fn is_empty(&self) -> bool {
        self.trim().is_empty()
    }

    fn extend(&mut self, other: Self) -> () {
        self.push_str(&other)
    }
}
