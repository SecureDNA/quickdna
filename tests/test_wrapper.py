from quickdna import DnaSequence, ProteinSequence

def test_translate():
    assert DnaSequence("AAAGGGAAA").translate(table=1) == ProteinSequence("KGK")

def test_equality():
    d1 = DnaSequence("aaa")
    d2 = DnaSequence("aaa")
    p1 = ProteinSequence("aaa")

    assert d1 == d2
    assert d1 != p1
    assert d2 != p1

def test_hash():
    d1 = DnaSequence("aaa")
    d2 = DnaSequence("aaa")
    p1 = ProteinSequence("aaa")

    assert hash(d1) == hash(d2)
    assert hash(d1) != hash(p1)

def test_indexing_slicing():
    d = DnaSequence("aaa")
    assert d[1] == ord("a")
    assert d[:] == DnaSequence("aaa")
    assert d[1:] == DnaSequence("aa")