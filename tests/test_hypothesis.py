from Bio.Seq import Seq
from hypothesis import given, strategies as st

from quickdna import DnaSequence

st_dna = st.text(alphabet=("a", "A", "t", "T", "c", "C", "g", "G"))

@given(st_dna)
def test_translate(dna):
    for table in range(1, 34):
        if table not in (7, 8, 17, 18, 19, 20):
            quickdna_translation = DnaSequence(dna).translate(table).seq
            biopython_translation = bytes(Seq(dna).translate(table=table))
            assert quickdna_translation == biopython_translation

@given(st_dna)
def test_reverse_complement(dna):
    quickdna_rc = DnaSequence(dna).reverse_complement().seq

    biopython_rc = Seq(dna).reverse_complement()
    biopython_rc = str(biopython_rc).upper().encode("ascii")

    assert quickdna_rc == biopython_rc
