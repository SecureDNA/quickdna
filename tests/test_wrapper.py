import pytest

from quickdna import DnaSequence, ProteinSequence


def test_translate():
    assert DnaSequence("AAAGGGAAA").translate(table=1) == ProteinSequence("KGK")


def test_translate_self():
    assert DnaSequence("AAAGGGAAA").translate_self_frames(table=1) == [
        ProteinSequence("KGK"),
        ProteinSequence("KG"),
        ProteinSequence("RE"),
    ]


def test_short_translate_self():
    assert DnaSequence("GGGG").translate_self_frames(table=1) == [
        ProteinSequence("G"),
        ProteinSequence("G"),
    ]

    assert DnaSequence("GGG").translate_self_frames(table=1) == [
        ProteinSequence("G"),
    ]

    assert DnaSequence("GG").translate_self_frames(table=1) == []
    assert DnaSequence("G").translate_self_frames(table=1) == []
    assert DnaSequence("").translate_self_frames(table=1) == []


def test_translate_all():
    assert DnaSequence("AAAGGGAAA").translate_all_frames(table=1) == [
        ProteinSequence("KGK"),
        ProteinSequence("KG"),
        ProteinSequence("RE"),
        ProteinSequence("FPF"),
        ProteinSequence("FP"),
        ProteinSequence("SL"),
    ]


def test_short_translate_all():
    assert DnaSequence("GGGG").translate_all_frames(table=1) == [
        ProteinSequence("G"),
        ProteinSequence("G"),
        ProteinSequence("P"),
        ProteinSequence("P"),
    ]

    assert DnaSequence("GGG").translate_all_frames(table=1) == [
        ProteinSequence("G"),
        ProteinSequence("P"),
    ]

    assert DnaSequence("GG").translate_all_frames(table=1) == []
    assert DnaSequence("G").translate_all_frames(table=1) == []
    assert DnaSequence("").translate_all_frames(table=1) == []


def test_equality():
    d1 = DnaSequence("aaa")
    d2 = DnaSequence("aaa")
    p1 = ProteinSequence("aaa")
    p2 = ProteinSequence("aaa")

    assert d1 == d2
    assert d1 != p1
    assert d2 != p1
    assert p1 == p2


def test_hash():
    d1 = DnaSequence("aaa")
    d2 = DnaSequence("aaa")
    p1 = ProteinSequence("aaa")
    p2 = ProteinSequence("aaa")

    assert hash(d1) == hash(d2)
    assert hash(d1) != hash(p1)
    assert hash(p1) == hash(p2)


def test_indexing_slicing_dna():
    d = DnaSequence("aaa")
    assert d[1] == ord("a")
    assert d[:] == DnaSequence("aaa")
    assert d[1:] == DnaSequence("aa")
    assert d[1:] != ProteinSequence("aa")


def test_indexing_slicing_protein():
    p = ProteinSequence("aaa")
    assert p[1] == ord("a")
    assert p[:] == ProteinSequence("aaa")
    assert p[1:] == ProteinSequence("aa")
    assert p[1:] != DnaSequence("aa")


def test_len_dna():
    d = DnaSequence("a" * 20)
    assert len(d) == 20


def test_len_protein():
    p = ProteinSequence("a" * 20)
    assert len(p) == 20


def test_add():
    d1 = DnaSequence("at")
    d2 = DnaSequence("cg")
    assert d1 + d2 == DnaSequence("atcg")

    p1 = ProteinSequence("KG")
    p2 = ProteinSequence("GK")
    assert p1 + p2 == ProteinSequence("KGGK")

    with pytest.raises(TypeError):
        d1 + p1  # type: ignore


def test_mul():
    d = DnaSequence("a") * 20
    assert d == DnaSequence("a" * 20)

    p = ProteinSequence("k") * 20
    assert p == ProteinSequence("k" * 20)


def test_iter():
    dna_src = "atcg" * 20
    assert list(DnaSequence(dna_src)) == list(dna_src)

    p_src = "kgkg" * 20
    assert list(ProteinSequence(p_src)) == list(p_src)
