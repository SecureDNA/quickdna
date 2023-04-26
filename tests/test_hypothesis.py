# Copyright 2021-2023 SecureDNA Stiftung (SecureDNA Foundation) <license@securedna.org>
# SPDX-License-Identifier: MIT OR Apache-2.0

"""
These are conformance tests that ensure the output of quickdna matches biopython
for randomly-generated strings of DNA.

Hypothesis is a property testing library, which means it takes some schema, in
our case valid DNA strings, and generates random data matching that schema in
order to test certain assertions hold.
"""

from Bio.Seq import Seq
from hypothesis import given, strategies as st

from quickdna import DnaSequence
from .utils import valid_tables

st_dna = st.text(alphabet=("a", "A", "t", "T", "c", "C", "g", "G", "n", "N"))


@given(st_dna)
def test_translate(dna):
    for table in valid_tables():
        quickdna_translation = DnaSequence(dna).translate(table).seq
        biopython_translation = bytes(Seq(dna).translate(table=table))
        assert quickdna_translation == biopython_translation


@given(st_dna)
def test_reverse_complement(dna):
    quickdna_rc = DnaSequence(dna).reverse_complement().seq

    biopython_rc = Seq(dna).reverse_complement()
    biopython_rc = str(biopython_rc).upper().encode("ascii")

    assert quickdna_rc == biopython_rc
