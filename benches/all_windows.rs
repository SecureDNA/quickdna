use criterion::{black_box, criterion_group, criterion_main, Criterion};
use rand::{rngs::OsRng, seq::SliceRandom};
use smallvec::SmallVec;

use quickdna::{
    BaseSequence, DnaSequence, Nucleotide, NucleotideLike, Nucleotides, TranslationTable,
};

static PROTEIN_WINDOW_LEN: usize = 20;
static DNA_WINDOW_LEN: usize = 42;

// Adapted from synthclient
#[derive(Debug, Clone)]
pub struct SequenceWindows {
    pub dna_sequence_length: usize,
    pub dna: Vec<String>,
    pub dna_reverse: Vec<String>,
    pub aa: Vec<Vec<String>>,
}

impl SequenceWindows {
    /// returns the combined length of all types of windows
    fn len(&self) -> usize {
        self.dna.len()
            + self.dna_reverse.len()
            + self
                .aa
                .iter()
                .map(|ele| ele.len())
                .reduce(|a, b| a + b)
                .unwrap_or(0)
    }

    /// Returns an iterator over the windows and their original position in the full FASTA
    /// The indexes returned by enumerate are not strictly monotonically increasing
    fn enumerate(&self) -> impl Iterator<Item = (usize, &str)> {
        self.dna
            .iter()
            .enumerate()
            .map(|(i, s)| (i, s.as_str()))
            .chain(
                self.dna_reverse
                    .iter()
                    .enumerate()
                    .map(|(index, value)| (self.dna_reverse.len() - index - 1, value.as_str())),
            )
            .chain(
                self.aa
                    .iter()
                    .enumerate()
                    .flat_map(move |(operation, frame_windows)| {
                        // this code is the reverse of QuickDNA translate_all_frames
                        // https://github.com/SecureDNA/quickdna/blob/main/src/rust_api.rs#L196
                        let dna_per_aa = 3;
                        frame_windows
                            .iter()
                            .enumerate()
                            .map(move |(index_within_frame, window)| {
                                if operation == 0 || operation == 1 || operation == 2 {
                                    // forward
                                    (dna_per_aa * index_within_frame + operation, window.as_str())
                                } else {
                                    // backwards
                                    let index = self.dna_sequence_length
                                        - ((index_within_frame + PROTEIN_WINDOW_LEN) * dna_per_aa
                                            + (operation % 3));
                                    (index, window.as_str())
                                }
                            })
                    }),
            )
    }

    /// Construct all windows from a DNA sequence
    fn from_dna(dna: &DnaSequence<Nucleotide>) -> SequenceWindows {
        let translations = dna.translate_all_frames(TranslationTable::Ncbi1);

        let windows = SequenceWindows {
            dna_sequence_length: dna.len(),
            dna: dna.windows(DNA_WINDOW_LEN).map(|s| s.to_string()).collect(),
            dna_reverse: dna
                .reverse_complement()
                .windows(DNA_WINDOW_LEN)
                .map(|s| s.to_string())
                .collect(),
            aa: translations
                .iter()
                .map(|t| {
                    t.windows(PROTEIN_WINDOW_LEN)
                        .map(|s| s.to_string())
                        .collect()
                })
                .collect(),
        };

        // amino acids
        windows
    }
}

#[derive(Debug, Clone)]
pub struct IterBasedSequenceWindows {
    dna_window_len: usize,
    protein_window_len: usize,
    dna: String,
    dna_rc: String,
    aas: SmallVec<[String; 6]>,
}

impl IterBasedSequenceWindows {
    fn from_dna(dna: &[Nucleotide], dna_window_len: usize, protein_window_len: usize) -> Self {
        let mut aas = SmallVec::new();
        let ncbi1 = TranslationTable::Ncbi1.to_fn();
        aas.extend(dna.self_reading_frames().into_iter().map(|codons| {
            let translated: Vec<_> = codons.map(ncbi1).collect();
            String::from_utf8(translated).unwrap()
        }));
        aas.extend(
            dna.reverse_complement()
                .self_reading_frames()
                .into_iter()
                .map(|codons| {
                    let translated: Vec<_> = codons.map(ncbi1).collect();
                    String::from_utf8(translated).unwrap()
                }),
        );

        let dna_rc = (&dna).reverse_complement().map(|n| n.to_ascii()).collect();
        let dna_rc = String::from_utf8(dna_rc).unwrap();

        let dna = dna.iter().map(|n| n.to_ascii()).collect();
        let dna = String::from_utf8(dna).unwrap();

        Self {
            dna_window_len,
            protein_window_len,
            dna,
            dna_rc,
            aas,
        }
    }

    fn enumerate(&self) -> impl Iterator<Item = (usize, &str)> {
        self.dna_windows()
            .chain(self.dna_rc_windows())
            .chain(self.aa_windows().flatten())
    }

    fn len(&self) -> usize {
        let aas_len: usize = self.aa_windows().map(|iter| iter.len()).sum();
        self.dna_windows().len() + self.dna_rc_windows().len() + aas_len
    }

    fn dna_windows(&self) -> impl ExactSizeIterator<Item = (usize, &str)> {
        Self::ascii_str_windows(&self.dna, self.dna_window_len).enumerate()
    }

    fn dna_rc_windows(&self) -> impl ExactSizeIterator<Item = (usize, &str)> {
        let first_dna_rc_idx = Self::num_windows(self.dna_rc.len(), self.dna_window_len) - 1;
        Self::ascii_str_windows(&self.dna_rc, self.dna_window_len)
            .enumerate()
            .map(move |(i, w)| (first_dna_rc_idx - i, w))
    }

    fn aa_windows(
        &self,
    ) -> impl Iterator<Item = impl ExactSizeIterator<Item = (usize, &str)> + '_> + '_ {
        // Thankfully, we don't actually need to worry about the edge-case where there are
        // fewer than 6 reading frames. That can only happen if the dna sequence is <5 bases
        // long, in which case the aas would be less than the protein window length and
        // therefore lack any windows, making the exact logic of window generation moot.
        (0..self.aas.len()).map(|frame_idx| {
            let indices = self.aa_window_indices(frame_idx);
            let windows = Self::ascii_str_windows(&self.aas[frame_idx], self.protein_window_len);
            indices.zip(windows)
        })
    }

    // Note: Imitates original logic so I can verify the iter-based window code is correct.
    fn aa_window_indices(&self, frame_idx: usize) -> impl ExactSizeIterator<Item = usize> {
        let dna_len = self.dna.len();
        let protein_window_len = self.protein_window_len;
        let num_windows = Self::num_windows(dna_len, protein_window_len);
        let is_reverse = frame_idx >= 3;
        let frame_offset = frame_idx % 3;
        (0..num_windows).map(move |i| {
            if is_reverse {
                dna_len - ((i + protein_window_len) * 3 + frame_offset)
            } else {
                3 * i + frame_offset
            }
        })
    }

    // Unfortunate, but needed as long as we're using strings, AFAICT
    fn ascii_str_windows(s: &str, window_len: usize) -> impl ExactSizeIterator<Item = &str> {
        let num_windows = Self::num_windows(s.len(), window_len);
        (0..num_windows).map(move |i| &s[i..i + window_len])
    }

    fn num_windows(data_len: usize, window_len: usize) -> usize {
        data_len + 1 - window_len.min(data_len + 1)
    }
}

pub fn criterion_benchmark(c: &mut Criterion) {
    let nucleotides = [Nucleotide::A, Nucleotide::T, Nucleotide::C, Nucleotide::G];
    let dna: Vec<_> = (0..99999)
        .map(|_| *nucleotides.choose(&mut OsRng).unwrap())
        .collect();
    let dna = DnaSequence::new(dna);

    // Sanity check that the iterator-based code behaves the same.
    let baseline_windows = SequenceWindows::from_dna(&dna);
    let new_windows =
        IterBasedSequenceWindows::from_dna(dna.as_slice(), DNA_WINDOW_LEN, PROTEIN_WINDOW_LEN);
    assert_eq!(new_windows.len(), baseline_windows.len());
    assert!(new_windows.enumerate().eq(baseline_windows.enumerate()));

    c.bench_function("DNASequence-based windows", |b| {
        b.iter(|| {
            let windows = SequenceWindows::from_dna(black_box(&dna));
            for window in windows.enumerate() {
                black_box(window);
            }
        })
    });

    c.bench_function("Iter-based windows", |b| {
        b.iter(|| {
            let windows = IterBasedSequenceWindows::from_dna(
                black_box(dna.as_slice()),
                DNA_WINDOW_LEN,
                PROTEIN_WINDOW_LEN,
            );
            for window in windows.enumerate() {
                black_box(window);
            }
        })
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
