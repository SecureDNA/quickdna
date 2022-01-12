# quickdna

[![PyPI](https://img.shields.io/pypi/v/quickdna?style=flat-square)](https://pypi.org/project/quickdna/)

quickdna is a simple, fast library for working with DNA sequences.

```python
# These are the two main library types. Unlike Biopython, DnaSequence and ProteinSequence are distinct,
# though they share a common BaseSequence base class
>>> from quickdna import DnaSequence, ProteinSequence

# Sequences can be constructed from strs or bytes, and are stored internally as ascii-encoded bytes
>>> d = DnaSequence("taatcaagactattcaaccaa")

# If no table=... argument is given, NBCI table 1 will be used by default...
>>> d.translate()
ProteinSequence(seq='*SRLFNQ')

# But any of the NCBI tables can be specified. A ValueError will be thrown for an invalid table.
>>> d.translate(table=22)
ProteinSequence(seq='**RLFNQ')

# This method will return a tuple of all possible reading frames (seq[:], seq[1:], and seq[2:])
>>> d.translate_all_frames()
(ProteinSequence(seq='*SRLFNQ'), ProteinSequence(seq='NQDYST'), ProteinSequence(seq='IKTIQP'))

# Sequences can be sliced just like regular strings, and return new sequence instances.
>>> d[3:9].translate()
ProteinSequence(seq='SR')

# This exists too!
>>> d[3:9].reverse_complement()
DnaSequence(seq='TCTTGA')
```

For regular DNA translation tasks, quickdna is faster than Biopython. (See `benchmarks/bench.py` for source)

task                                       | time             | comparison
-------------------------------------------|------------------|-----------
translate_quickdna(small_genome)           | 0.00306ms / iter |
translate_biopython(small_genome)          | 0.05834ms / iter | 1908.90%
translate_quickdna(covid_genome)           | 0.02959ms / iter |
translate_biopython(covid_genome)          | 3.54413ms / iter | 11979.10%
reverse_complement_quickdna(small_genome)  | 0.00238ms / iter |
reverse_complement_biopython(small_genome) | 0.00398ms / iter | 167.24%
reverse_complement_quickdna(covid_genome)  | 0.02409ms / iter |
reverse_complement_biopython(covid_genome) | 0.02928ms / iter | 121.55%

