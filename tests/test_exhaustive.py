"""
This test simply ensures every possible codon & table combination matches Biopython
"""

from Bio.Seq import Seq

from quickdna import DnaSequence
from .utils import valid_tables

NUCLEOTIDES = "ATCGN"


def test_exhaustive():
    for table in valid_tables():
        for a in NUCLEOTIDES:
            for b in NUCLEOTIDES:
                for c in NUCLEOTIDES:
                    dna = f"{a}{b}{c}"
                    qd_trans = chr(DnaSequence(dna).translate(table=table)[0])
                    bp_trans = str(Seq(dna).translate(table=table))[0]

                    assert qd_trans == bp_trans, dna
