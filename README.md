# quickdna

quickdna is a simple, fast library for working with DNA sequences.

```python
>>> from quickdna import DnaSequence, ProteinSequence
>>> d = DnaSequence("taatcaagactattcaaccaa")
>>> d.translate()
ProteinSequence(seq='*SRLFNQ')
>>> d.translate(table=22)
ProteinSequence(seq='**RLFNQ')
>>> d.translate_all_frames()
(ProteinSequence(seq='*SRLFNQ'), ProteinSequence(seq='NQDYST'), ProteinSequence(seq='IKTIQP'))
>>> d[3:9].translate()
ProteinSequence(seq='SR')
>>> d[3:9].reverse_complement()
DnaSequence(seq='TCTTGA')
```

quickdna is much faster than Biopython for regular DNA translation tasks.

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

