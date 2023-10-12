use std::cmp::Ordering;

use crate::{Nucleotide, NucleotideLike};

pub struct ForwardCanonicalized<I> {
    permutation: u16,
    inner: I,
}

impl<N, I> Iterator for ForwardCanonicalized<I>
where
    N: NucleotideLike,
    I: Iterator<Item = N>,
{
    type Item = Nucleotide;

    fn next(&mut self) -> Option<Self::Item> {
        let nuc = self.inner.next()?.bits() as u16;

        for (i, out_nuc) in Nucleotide::ALL.into_iter().enumerate() {
            let mask = (nuc & 0x000F) << (4 * i);
            if self.permutation & mask != 0 {
                let nonmask_bits = !(0x000F << (4 * i));
                self.permutation &= mask | nonmask_bits;
                return Some(out_nuc);
            }
        }
        panic!("Unknown bug: Overconstrained permutation. Maybe the nucleotide's binary representation is invalid?");
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.inner.size_hint()
    }
}

impl<N: NucleotideLike, I: ExactSizeIterator<Item = N>> ExactSizeIterator
    for ForwardCanonicalized<I>
{
}

// Canonicalizes only in one direction
// Note: This is only guaranteed to work with Nucleotide and NucleotideAmbiguous iters.
pub fn forward_canonicalize<N, I>(iterable: I) -> ForwardCanonicalized<I::IntoIter>
where
    N: NucleotideLike,
    I: IntoIterator<Item = N>,
{
    ForwardCanonicalized {
        inner: iterable.into_iter(),
        permutation: !0,
    }
}

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

pub struct Canonicalized<I>(
    LexicalMin<ForwardCanonicalized<I>, ForwardCanonicalized<std::iter::Rev<I>>>,
);

impl<I> Iterator for Canonicalized<I>
where
    I: DoubleEndedIterator,
    I::Item: NucleotideLike,
{
    type Item = Nucleotide;

    fn next(&mut self) -> Option<Self::Item> {
        self.0.next()
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.0.size_hint()
    }
}

impl<I> ExactSizeIterator for Canonicalized<I>
where
    I: ExactSizeIterator + DoubleEndedIterator,
    I::Item: NucleotideLike,
{
}

pub fn canonicalize<N, I>(iterable: I) -> Canonicalized<I::IntoIter>
where
    N: NucleotideLike,
    I: IntoIterator<Item = N>,
    I::IntoIter: DoubleEndedIterator + Clone,
{
    let iter = iterable.into_iter();
    let forward_canonical = forward_canonicalize(iter.clone());
    let reverse_canonical = forward_canonicalize(iter.rev());
    Canonicalized(LexicalMin::new(forward_canonical, reverse_canonical))
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct Remap(pub [Nucleotide; 4]);

impl Remap {
    pub fn apply(self, nuc: Nucleotide) -> Nucleotide {
        self.0[(nuc as u8).ilog2() as usize]
    }

    pub fn remap(self, dna: &[Nucleotide]) -> Vec<Nucleotide> {
        dna.iter().copied().map(|nuc| self.apply(nuc)).collect()
    }
}

pub const PERMUTATIONS: [Remap; 24] = {
    use Nucleotide::{A, C, G, T};
    use Remap as R;
    [
        R([A, T, C, G]),
        R([A, T, G, C]),
        R([A, C, T, G]),
        R([A, C, G, T]),
        R([A, G, T, C]),
        R([A, G, C, T]),
        R([T, A, C, G]),
        R([T, A, G, C]),
        R([T, C, A, G]),
        R([T, C, G, A]),
        R([T, G, A, C]),
        R([T, G, C, A]),
        R([C, A, T, G]),
        R([C, A, G, T]),
        R([C, T, A, G]),
        R([C, T, G, A]),
        R([C, G, A, T]),
        R([C, G, T, A]),
        R([G, A, T, C]),
        R([G, A, C, T]),
        R([G, T, A, C]),
        R([G, T, C, A]),
        R([G, C, A, T]),
        R([G, C, T, A]),
    ]
};

#[cfg(test)]
mod test {
    use super::*;

    use quickcheck::{quickcheck, Arbitrary, Gen};

    use crate::{BaseSequence, DnaSequenceAmbiguous, DnaSequenceStrict, NucleotideAmbiguous};

    fn dna(dna: &str) -> DnaSequenceStrict {
        dna.parse().unwrap()
    }

    fn strict_canon(src_dna: &str) -> String {
        let src_dna: DnaSequenceStrict = src_dna.parse().unwrap();
        let canonical = canonicalize(src_dna.as_slice().iter().copied()).collect();
        DnaSequenceStrict::new(canonical).to_string()
    }

    fn canon(src_dna: &str) -> String {
        let src_dna: DnaSequenceAmbiguous = src_dna.parse().unwrap();
        let canonical = canonicalize(src_dna.as_slice().iter().copied()).collect();
        DnaSequenceStrict::new(canonical).to_string()
    }

    #[test]
    fn sanity_check_strict_canonicalize() {
        // Testing forward vs reverse canonicalization
        assert_eq!(canon("CCA"), "AAT");
        assert_eq!(canon("ACC"), "AAT");
        assert_eq!(canon("TTGT"), "AATA");
        assert_eq!(canon("TGTT"), "AATA");
        // Anything unambiguous that begins with AATCG and ends with two different nucleotides
        // is already canonical
        assert_eq!(
            canon("AATCGATCGTAGCATGACTGACTGCTAATGTCGTAGCGTAGT"),
            "AATCGATCGTAGCATGACTGACTGCTAATGTCGTAGCGTAGT"
        );
    }

    #[test]
    fn sanity_check_ambiguous_canonicalize() {
        // Testing forward vs reverse canonicalization
        assert_eq!(canon("CCA"), "AAT");
        assert_eq!(canon("ACC"), "AAT");
        assert_eq!(canon("TTGT"), "AATA");
        assert_eq!(canon("TGTT"), "AATA");
        // Testing ambiguity codes; Y=TC, M=AC
        assert_eq!(canon("YAM"), "AAT"); // reversed is earlier
        assert_eq!(canon("YTM"), "AAT");
        assert_eq!(canon("YCM"), "AAA");
        assert_eq!(canon("YGM"), "ATA");
        // Anything unambiguous that begins with AATCG and ends with two different nucleotides
        // is already canonical
        assert_eq!(
            canon("AATCGATCGTAGCATGACTGACTGCTAATGTCGTAGCGTAGT"),
            "AATCGATCGTAGCATGACTGACTGCTAATGTCGTAGCGTAGT"
        );
        // Something with a little bit of everyhing; W=AT, M=AC, N=ATCG
        assert_eq!(
            canon("ACTGAWMTGCTGACGAGTANGTCGAGGT"),
            "ATTCTGATACATCTGCTAGTACACTAGC"
        );
    }

    #[test]
    fn exhaustively_check_small_strict_dna() {
        assert_eq!(strict_canon(""), "");

        assert_eq!(strict_canon("A"), "A");
        assert_eq!(strict_canon("T"), "A");
        assert_eq!(strict_canon("C"), "A");
        assert_eq!(strict_canon("G"), "A");

        assert_eq!(strict_canon("AA"), "AA");
        assert_eq!(strict_canon("AT"), "AT");
        assert_eq!(strict_canon("AC"), "AT");
        assert_eq!(strict_canon("AG"), "AT");
        assert_eq!(strict_canon("TA"), "AT");
        assert_eq!(strict_canon("TT"), "AA");
        assert_eq!(strict_canon("TC"), "AT");
        assert_eq!(strict_canon("TG"), "AT");
        assert_eq!(strict_canon("CA"), "AT");
        assert_eq!(strict_canon("CT"), "AT");
        assert_eq!(strict_canon("CC"), "AA");
        assert_eq!(strict_canon("CG"), "AT");
        assert_eq!(strict_canon("GA"), "AT");
        assert_eq!(strict_canon("GT"), "AT");
        assert_eq!(strict_canon("GC"), "AT");
        assert_eq!(strict_canon("GG"), "AA");

        assert_eq!(strict_canon("AAA"), "AAA");
        assert_eq!(strict_canon("AAT"), "AAT");
        assert_eq!(strict_canon("AAC"), "AAT");
        assert_eq!(strict_canon("AAG"), "AAT");
        assert_eq!(strict_canon("ATA"), "ATA");
        assert_eq!(strict_canon("ATT"), "AAT");
        assert_eq!(strict_canon("ATC"), "ATC");
        assert_eq!(strict_canon("ATG"), "ATC");
        assert_eq!(strict_canon("ACA"), "ATA");
        assert_eq!(strict_canon("ACT"), "ATC");
        assert_eq!(strict_canon("ACC"), "AAT");
        assert_eq!(strict_canon("ACG"), "ATC");
        assert_eq!(strict_canon("AGA"), "ATA");
        assert_eq!(strict_canon("AGT"), "ATC");
        assert_eq!(strict_canon("AGC"), "ATC");
        assert_eq!(strict_canon("AGG"), "AAT");

        assert_eq!(strict_canon("TAA"), "AAT");
        assert_eq!(strict_canon("TAT"), "ATA");
        assert_eq!(strict_canon("TAC"), "ATC");
        assert_eq!(strict_canon("TAG"), "ATC");
        assert_eq!(strict_canon("TTA"), "AAT");
        assert_eq!(strict_canon("TTT"), "AAA");
        assert_eq!(strict_canon("TTC"), "AAT");
        assert_eq!(strict_canon("TTG"), "AAT");
        assert_eq!(strict_canon("TCA"), "ATC");
        assert_eq!(strict_canon("TCT"), "ATA");
        assert_eq!(strict_canon("TCC"), "AAT");
        assert_eq!(strict_canon("TCG"), "ATC");
        assert_eq!(strict_canon("TGA"), "ATC");
        assert_eq!(strict_canon("TGT"), "ATA");
        assert_eq!(strict_canon("TGC"), "ATC");
        assert_eq!(strict_canon("TGG"), "AAT");

        assert_eq!(strict_canon("CAA"), "AAT");
        assert_eq!(strict_canon("CAT"), "ATC");
        assert_eq!(strict_canon("CAC"), "ATA");
        assert_eq!(strict_canon("CAG"), "ATC");
        assert_eq!(strict_canon("CTA"), "ATC");
        assert_eq!(strict_canon("CTT"), "AAT");
        assert_eq!(strict_canon("CTC"), "ATA");
        assert_eq!(strict_canon("CTG"), "ATC");
        assert_eq!(strict_canon("CCA"), "AAT");
        assert_eq!(strict_canon("CCT"), "AAT");
        assert_eq!(strict_canon("CCC"), "AAA");
        assert_eq!(strict_canon("CCG"), "AAT");
        assert_eq!(strict_canon("CGA"), "ATC");
        assert_eq!(strict_canon("CGT"), "ATC");
        assert_eq!(strict_canon("CGC"), "ATA");
        assert_eq!(strict_canon("CGG"), "AAT");

        assert_eq!(strict_canon("GAA"), "AAT");
        assert_eq!(strict_canon("GAT"), "ATC");
        assert_eq!(strict_canon("GAC"), "ATC");
        assert_eq!(strict_canon("GAG"), "ATA");
        assert_eq!(strict_canon("GTA"), "ATC");
        assert_eq!(strict_canon("GTT"), "AAT");
        assert_eq!(strict_canon("GTC"), "ATC");
        assert_eq!(strict_canon("GTG"), "ATA");
        assert_eq!(strict_canon("GCA"), "ATC");
        assert_eq!(strict_canon("GCT"), "ATC");
        assert_eq!(strict_canon("GCC"), "AAT");
        assert_eq!(strict_canon("GCG"), "ATA");
        assert_eq!(strict_canon("GGA"), "AAT");
        assert_eq!(strict_canon("GGT"), "AAT");
        assert_eq!(strict_canon("GGC"), "AAT");
        assert_eq!(strict_canon("GGG"), "AAA");
    }

    #[test]
    fn exhaustively_check_small_ambiguous_dna() {
        assert_eq!(canon(""), "");

        assert_eq!(canon("A"), "A");
        assert_eq!(canon("T"), "A");
        assert_eq!(canon("C"), "A");
        assert_eq!(canon("G"), "A");
        assert_eq!(canon("W"), "A");
        assert_eq!(canon("M"), "A");
        assert_eq!(canon("R"), "A");
        assert_eq!(canon("Y"), "A");
        assert_eq!(canon("S"), "A");
        assert_eq!(canon("K"), "A");
        assert_eq!(canon("B"), "A");
        assert_eq!(canon("V"), "A");
        assert_eq!(canon("D"), "A");
        assert_eq!(canon("H"), "A");
        assert_eq!(canon("N"), "A");

        assert_eq!(canon("AA"), "AA");
        assert_eq!(canon("AT"), "AT");
        assert_eq!(canon("AC"), "AT");
        assert_eq!(canon("AG"), "AT");
        assert_eq!(canon("AW"), "AA");
        assert_eq!(canon("AM"), "AA");
        assert_eq!(canon("AR"), "AA");
        assert_eq!(canon("AY"), "AT");
        assert_eq!(canon("AS"), "AT");
        assert_eq!(canon("AK"), "AT");
        assert_eq!(canon("AB"), "AT");
        assert_eq!(canon("AV"), "AA");
        assert_eq!(canon("AD"), "AA");
        assert_eq!(canon("AH"), "AA");
        assert_eq!(canon("AN"), "AA");

        assert_eq!(canon("TA"), "AT");
        assert_eq!(canon("TT"), "AA");
        assert_eq!(canon("TC"), "AT");
        assert_eq!(canon("TG"), "AT");
        assert_eq!(canon("TW"), "AA");
        assert_eq!(canon("TM"), "AT");
        assert_eq!(canon("TR"), "AT");
        assert_eq!(canon("TY"), "AA");
        assert_eq!(canon("TS"), "AT");
        assert_eq!(canon("TK"), "AA");
        assert_eq!(canon("TB"), "AA");
        assert_eq!(canon("TV"), "AT");
        assert_eq!(canon("TD"), "AA");
        assert_eq!(canon("TH"), "AA");
        assert_eq!(canon("TN"), "AA");

        assert_eq!(canon("CA"), "AT");
        assert_eq!(canon("CT"), "AT");
        assert_eq!(canon("CC"), "AA");
        assert_eq!(canon("CG"), "AT");
        assert_eq!(canon("CW"), "AT");
        assert_eq!(canon("CM"), "AA");
        assert_eq!(canon("CR"), "AT");
        assert_eq!(canon("CY"), "AA");
        assert_eq!(canon("CS"), "AA");
        assert_eq!(canon("CK"), "AT");
        assert_eq!(canon("CB"), "AA");
        assert_eq!(canon("CV"), "AA");
        assert_eq!(canon("CD"), "AT");
        assert_eq!(canon("CH"), "AA");
        assert_eq!(canon("CN"), "AA");

        assert_eq!(canon("GA"), "AT");
        assert_eq!(canon("GT"), "AT");
        assert_eq!(canon("GC"), "AT");
        assert_eq!(canon("GG"), "AA");
        assert_eq!(canon("GW"), "AT");
        assert_eq!(canon("GM"), "AT");
        assert_eq!(canon("GR"), "AA");
        assert_eq!(canon("GY"), "AT");
        assert_eq!(canon("GS"), "AA");
        assert_eq!(canon("GK"), "AA");
        assert_eq!(canon("GB"), "AA");
        assert_eq!(canon("GV"), "AA");
        assert_eq!(canon("GD"), "AA");
        assert_eq!(canon("GH"), "AT");
        assert_eq!(canon("GN"), "AA");

        assert_eq!(canon("WA"), "AA");
        assert_eq!(canon("WT"), "AA");
        assert_eq!(canon("WC"), "AT");
        assert_eq!(canon("WG"), "AT");
        assert_eq!(canon("WW"), "AA");
        assert_eq!(canon("WM"), "AA");
        assert_eq!(canon("WR"), "AA");
        assert_eq!(canon("WY"), "AA");
        assert_eq!(canon("WS"), "AT");
        assert_eq!(canon("WK"), "AA");
        assert_eq!(canon("WB"), "AA");
        assert_eq!(canon("WV"), "AA");
        assert_eq!(canon("WD"), "AA");
        assert_eq!(canon("WH"), "AA");
        assert_eq!(canon("WN"), "AA");

        assert_eq!(canon("MA"), "AA");
        assert_eq!(canon("MT"), "AT");
        assert_eq!(canon("MC"), "AA");
        assert_eq!(canon("MG"), "AT");
        assert_eq!(canon("MW"), "AA");
        assert_eq!(canon("MM"), "AA");
        assert_eq!(canon("MR"), "AA");
        assert_eq!(canon("MY"), "AA");
        assert_eq!(canon("MS"), "AA");
        assert_eq!(canon("MK"), "AT");
        assert_eq!(canon("MB"), "AA");
        assert_eq!(canon("MV"), "AA");
        assert_eq!(canon("MD"), "AA");
        assert_eq!(canon("MH"), "AA");
        assert_eq!(canon("MN"), "AA");

        assert_eq!(canon("RA"), "AA");
        assert_eq!(canon("RT"), "AT");
        assert_eq!(canon("RC"), "AT");
        assert_eq!(canon("RG"), "AA");
        assert_eq!(canon("RW"), "AA");
        assert_eq!(canon("RM"), "AA");
        assert_eq!(canon("RR"), "AA");
        assert_eq!(canon("RY"), "AT");
        assert_eq!(canon("RS"), "AA");
        assert_eq!(canon("RK"), "AA");
        assert_eq!(canon("RB"), "AA");
        assert_eq!(canon("RV"), "AA");
        assert_eq!(canon("RD"), "AA");
        assert_eq!(canon("RH"), "AA");
        assert_eq!(canon("RN"), "AA");

        assert_eq!(canon("YA"), "AT");
        assert_eq!(canon("YT"), "AA");
        assert_eq!(canon("YC"), "AA");
        assert_eq!(canon("YG"), "AT");
        assert_eq!(canon("YW"), "AA");
        assert_eq!(canon("YM"), "AA");
        assert_eq!(canon("YR"), "AT");
        assert_eq!(canon("YY"), "AA");
        assert_eq!(canon("YS"), "AA");
        assert_eq!(canon("YK"), "AA");
        assert_eq!(canon("YB"), "AA");
        assert_eq!(canon("YV"), "AA");
        assert_eq!(canon("YD"), "AA");
        assert_eq!(canon("YH"), "AA");
        assert_eq!(canon("YN"), "AA");

        assert_eq!(canon("SA"), "AT");
        assert_eq!(canon("ST"), "AT");
        assert_eq!(canon("SC"), "AA");
        assert_eq!(canon("SG"), "AA");
        assert_eq!(canon("SW"), "AT");
        assert_eq!(canon("SM"), "AA");
        assert_eq!(canon("SR"), "AA");
        assert_eq!(canon("SY"), "AA");
        assert_eq!(canon("SS"), "AA");
        assert_eq!(canon("SK"), "AA");
        assert_eq!(canon("SB"), "AA");
        assert_eq!(canon("SV"), "AA");
        assert_eq!(canon("SD"), "AA");
        assert_eq!(canon("SH"), "AA");
        assert_eq!(canon("SN"), "AA");

        assert_eq!(canon("KA"), "AT");
        assert_eq!(canon("KT"), "AA");
        assert_eq!(canon("KC"), "AT");
        assert_eq!(canon("KG"), "AA");
        assert_eq!(canon("KW"), "AA");
        assert_eq!(canon("KM"), "AT");
        assert_eq!(canon("KR"), "AA");
        assert_eq!(canon("KY"), "AA");
        assert_eq!(canon("KS"), "AA");
        assert_eq!(canon("KK"), "AA");
        assert_eq!(canon("KB"), "AA");
        assert_eq!(canon("KV"), "AA");
        assert_eq!(canon("KD"), "AA");
        assert_eq!(canon("KH"), "AA");
        assert_eq!(canon("KN"), "AA");

        assert_eq!(canon("BA"), "AT");
        assert_eq!(canon("BT"), "AA");
        assert_eq!(canon("BC"), "AA");
        assert_eq!(canon("BG"), "AA");
        assert_eq!(canon("BW"), "AA");
        assert_eq!(canon("BM"), "AA");
        assert_eq!(canon("BR"), "AA");
        assert_eq!(canon("BY"), "AA");
        assert_eq!(canon("BS"), "AA");
        assert_eq!(canon("BK"), "AA");
        assert_eq!(canon("BB"), "AA");
        assert_eq!(canon("BV"), "AA");
        assert_eq!(canon("BD"), "AA");
        assert_eq!(canon("BH"), "AA");
        assert_eq!(canon("BN"), "AA");

        assert_eq!(canon("VA"), "AA");
        assert_eq!(canon("VT"), "AT");
        assert_eq!(canon("VC"), "AA");
        assert_eq!(canon("VG"), "AA");
        assert_eq!(canon("VW"), "AA");
        assert_eq!(canon("VM"), "AA");
        assert_eq!(canon("VR"), "AA");
        assert_eq!(canon("VY"), "AA");
        assert_eq!(canon("VS"), "AA");
        assert_eq!(canon("VK"), "AA");
        assert_eq!(canon("VB"), "AA");
        assert_eq!(canon("VV"), "AA");
        assert_eq!(canon("VD"), "AA");
        assert_eq!(canon("VH"), "AA");
        assert_eq!(canon("VN"), "AA");

        assert_eq!(canon("DA"), "AA");
        assert_eq!(canon("DT"), "AA");
        assert_eq!(canon("DC"), "AT");
        assert_eq!(canon("DG"), "AA");
        assert_eq!(canon("DW"), "AA");
        assert_eq!(canon("DM"), "AA");
        assert_eq!(canon("DR"), "AA");
        assert_eq!(canon("DY"), "AA");
        assert_eq!(canon("DS"), "AA");
        assert_eq!(canon("DK"), "AA");
        assert_eq!(canon("DB"), "AA");
        assert_eq!(canon("DV"), "AA");
        assert_eq!(canon("DD"), "AA");
        assert_eq!(canon("DH"), "AA");
        assert_eq!(canon("DN"), "AA");

        assert_eq!(canon("HA"), "AA");
        assert_eq!(canon("HT"), "AA");
        assert_eq!(canon("HC"), "AA");
        assert_eq!(canon("HG"), "AT");
        assert_eq!(canon("HW"), "AA");
        assert_eq!(canon("HM"), "AA");
        assert_eq!(canon("HR"), "AA");
        assert_eq!(canon("HY"), "AA");
        assert_eq!(canon("HS"), "AA");
        assert_eq!(canon("HK"), "AA");
        assert_eq!(canon("HB"), "AA");
        assert_eq!(canon("HV"), "AA");
        assert_eq!(canon("HD"), "AA");
        assert_eq!(canon("HH"), "AA");
        assert_eq!(canon("HN"), "AA");

        assert_eq!(canon("NA"), "AA");
        assert_eq!(canon("NT"), "AA");
        assert_eq!(canon("NC"), "AA");
        assert_eq!(canon("NG"), "AA");
        assert_eq!(canon("NW"), "AA");
        assert_eq!(canon("NM"), "AA");
        assert_eq!(canon("NR"), "AA");
        assert_eq!(canon("NY"), "AA");
        assert_eq!(canon("NS"), "AA");
        assert_eq!(canon("NK"), "AA");
        assert_eq!(canon("NB"), "AA");
        assert_eq!(canon("NV"), "AA");
        assert_eq!(canon("ND"), "AA");
        assert_eq!(canon("NH"), "AA");
        assert_eq!(canon("NN"), "AA");
    }

    #[test]
    fn sanity_check_remap() {
        use Nucleotide::{A, C, G, T};
        let r = Remap([C, T, G, A]);
        let d = dna("ATCGAATTCCGG");
        let remapped = d.iter().map(|n| r.apply(n));
        assert!(remapped.eq(dna("CTGACCTTGGAA").iter()));
    }

    #[derive(Clone)]
    struct SemiAmbiguousDna(DnaSequenceAmbiguous);

    impl std::fmt::Debug for SemiAmbiguousDna {
        fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
            f.debug_tuple("SemiAmbiguousDna")
                .field(&self.0.to_string())
                .finish()
        }
    }

    impl Arbitrary for SemiAmbiguousDna {
        fn arbitrary(g: &mut Gen) -> Self {
            let dna: Vec<Nucleotide> = Arbitrary::arbitrary(g);
            let mut dna: Vec<_> = dna.into_iter().map(NucleotideAmbiguous::from).collect();
            // At 5 ambiguities, there may be as many as 1K expansions... I'm hesitant to go higher than that.
            let num_ambiguities = *g.choose(&[0, 1, 2, 3, 4, 5]).unwrap();
            if !dna.is_empty() {
                for _ in 0..num_ambiguities {
                    let i = usize::arbitrary(g) % dna.len();
                    dna[i] = NucleotideAmbiguous::arbitrary(g);
                }
            }
            Self(DnaSequenceAmbiguous::new(dna))
        }

        fn shrink(&self) -> Box<dyn Iterator<Item = Self>> {
            let base = self.0.as_slice().to_vec();
            let elements_removed = (0..base.len()).map(move |i| {
                let mut dna = base.clone();
                dna.swap_remove(i);
                Self(DnaSequenceAmbiguous::new(dna))
            });
            Box::new(elements_removed)
        }
    }

    fn strict_forward_canonicalize_reference_implementation(
        src_dna: &DnaSequenceStrict,
    ) -> Vec<Nucleotide> {
        PERMUTATIONS
            .iter()
            .map(|p| p.remap(src_dna.as_slice()))
            .min()
            .unwrap()
    }

    fn strict_canonicalize_reference_implementation(src_dna: DnaSequenceStrict) -> Vec<Nucleotide> {
        let fw_canonical = strict_forward_canonicalize_reference_implementation(&src_dna);
        let mut rev_dna = src_dna.as_slice().to_vec();
        rev_dna.reverse();
        let rev_dna = DnaSequenceStrict::new(rev_dna);
        let rev_canonical = strict_forward_canonicalize_reference_implementation(&rev_dna);
        fw_canonical.min(rev_canonical)
    }

    fn amb_forward_canonicalize_reference_implementation(
        src_dna: &DnaSequenceAmbiguous,
    ) -> Vec<Nucleotide> {
        src_dna
            .expansions()
            .flat_map(|expansion| PERMUTATIONS.iter().map(move |p| p.remap(&expansion)))
            .min()
            .unwrap()
    }

    fn amb_canonicalize_reference_implementation(src_dna: DnaSequenceAmbiguous) -> Vec<Nucleotide> {
        let fw_canonical = amb_forward_canonicalize_reference_implementation(&src_dna);
        let mut rev_dna = src_dna.as_slice().to_vec();
        rev_dna.reverse();
        let rev_dna = DnaSequenceAmbiguous::new(rev_dna);
        let rev_canonical = amb_forward_canonicalize_reference_implementation(&rev_dna);
        fw_canonical.min(rev_canonical)
    }

    quickcheck! {
        fn amb_canonicalize_matches_simple_definition(src_dna: SemiAmbiguousDna) -> bool {
            let actual: Vec<_> = canonicalize(src_dna.0.as_slice().iter().copied()).collect();
            let expected = amb_canonicalize_reference_implementation(src_dna.0);
            actual == expected
        }

        fn strict_canonicalize_matches_simple_definition(src_dna: DnaSequenceStrict) -> bool {
            let actual: Vec<_> = canonicalize(src_dna.as_slice().iter().copied()).collect();
            let expected = strict_canonicalize_reference_implementation(src_dna);
            actual == expected
        }
    }
}
