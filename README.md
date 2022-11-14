# quickdna

[![PyPI](https://img.shields.io/pypi/v/quickdna?style=flat-square)](https://pypi.org/project/quickdna/)

Quickdna is a simple, fast library for working with DNA sequences. It is up to 100x faster than Biopython for some
translation tasks, in part because it uses a native Rust module (via PyO3) for the translation. However, it exposes
an easy-to-use, type-annotated API that should still feel familiar for Biopython users.

```python
# These are the two main library types. Unlike Biopython, DnaSequence and
# ProteinSequence are distinct, though they share a common BaseSequence base class
>>> from quickdna import DnaSequence, ProteinSequence

# Sequences can be constructed from strs or bytes, and are stored internally as
# ascii-encoded bytes.
>>> d = DnaSequence("taatcaagactattcaaccaa")

# Sequences can be sliced just like regular strings, and return new sequence instances.
>>> d[3:9]
DnaSequence(seq='tcaaga')

# many other Python operations are supported on sequences as well: len, iter,
# ==, hash, concatenation with +, * a constant, etc. These operations are typed
# when appropriate and will not allow you to concatenate a ProteinSequence to a
# DnaSequence, for example

# DNA sequences can be easily translated to protein sequences with `translate()`.
# If no table=... argument is given, NBCI table 1 will be used by default...
>>> d.translate()
ProteinSequence(seq='*SRLFNQ')

# ...but any of the NCBI tables can be specified. A ValueError will be thrown
# for an invalid table.
>>> d.translate(table=22)
ProteinSequence(seq='**RLFNQ')

# This exists too! It's somewhat faster than Biopython, but not as dramatically as
# `translate()`
>>> d[3:9].reverse_complement()
DnaSequence(seq='TCTTGA')

# This method will return a list of all (up to 6) possible translated reading frames:
# (seq[:], seq[1:], seq[2:], seq.reverse_complement()[:], ...)
>>> d.translate_all_frames()
(ProteinSequence(seq='*SRLFNQ'), ProteinSequence(seq='NQDYST'),
ProteinSequence(seq='IKTIQP'), ProteinSequence(seq='LVE*S*L'),
ProteinSequence(seq='WLNSLD'), ProteinSequence(seq='G*IVLI'))

# translate_all_frames will return less than 6 frames for sequences of len < 5
>>> len(DnaSequence("AAAA").translate_all_frames())
4
>>> len(DnaSequence("AA").translate_all_frames())
0

# There is a similar method, `translate_self_frames`, that only returns the
# (up to 3) translated frames for this direction, without the reverse complement

# The IUPAC ambiguity codes are supported as well.
# Codons with N will translate to a specific amino acid if it is unambiguous,
# such as GGN -> G, or the ambiguous amino acid code 'X' if there are multiple
# possible translations.
>>> DnaSequence("GGNATN").translate()
ProteinSequence(seq='GX')

# The fine-grained ambiguity codes like "R = A or G" are accepted too, and
# translation results are the same as Biopython. In the output, amino acid
# ambiguity code 'B' means "either asparagine or aspartic acid" (N or D).
>>> DnaSequence("RAT").translate()
ProteinSequence(seq='B')

# To disallow ambiguity codes in translation, try: `.translate(strict=True)`
```

## Benchmarks

For regular DNA translation tasks, quickdna is faster than Biopython. (See `benchmarks/bench.py` for source).
Machines and workloads vary, however -- always benchmark!

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

## Should you use quickdna?

* Quickdna pros
  * It's quick!
  * It's simple and small.
  * It has type annotations, including a `py.typed` marker file for checkers like MyPy or VSCode's PyRight.
  * It makes a type distinction between DNA and protein sequences, preventing confusion.
* Quickdna cons:
  * It's newer and less battle-tested than Biopython.
  * It's not yet 1.0 -- the API is liable to change in the future.
  * It doesn't support reading FASTA files or many of the other tasks Biopython can do,
    so you'll probably end up still using Biopython or something else to do those tasks.

## Installation

Quickdna has prebuilt wheels for Linux (manylinux2010), OSX, and Windows available [on PyPi](https://pypi.org/project/quickdna/).

## Development

Quickdna uses `PyO3` and `maturin` to build and upload the wheels, and `poetry` for handling dependencies. This is handled via
a `Justfile`, which requires [Just](https://github.com/casey/just), a command-runner similar to `make`.

### Poetry

You can install poetry from https://python-poetry.org, and it will handle the other python dependencies.

### Just

You can install `Just` with `cargo install just`, and then run it in the project directory to get a list of commands.

### Flamegraphs

The `just profile` command requires [cargo-flamegraph](https://github.com/flamegraph-rs/flamegraph), please see that repository for installation instructions.
