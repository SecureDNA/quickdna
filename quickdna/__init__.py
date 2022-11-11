import typing as ty

from .quickdna import _translate, _translate_strict, _reverse_complement, _reverse_complement_strict  # type: ignore

T = ty.TypeVar("T", bound="BaseSequence")


def ensure_bytes(str_or_bytes: ty.Union[str, bytes]) -> bytes:
    """Encodes str as ascii (in strict mode), passing bytes along unchanged"""

    if isinstance(str_or_bytes, str):
        str_or_bytes = str_or_bytes.encode("ascii", "strict")
    return str_or_bytes


class BaseSequence:
    """Base class for DNA and Protein sequences."""

    def __init__(self, seq: ty.Union[str, bytes]) -> None:
        """
        Constructs a sequence from an input str or bytes.

        Raises UnicodeEncodeError if the input seq is a str with non-ascii characters.
        Does no other validation -- ValueError may be thrown later by `DnaSequence.translate()`,
        for example, if the sequence contains invalid characters like 'Q'.
        """

        self._seq = ensure_bytes(seq)

    @property
    def seq(self) -> bytes:
        return self._seq

    @ty.overload
    def __getitem__(self, __i: ty.SupportsIndex) -> int:
        ...

    @ty.overload
    def __getitem__(self: T, __s: slice) -> T:
        ...

    def __getitem__(self, key):
        if isinstance(key, ty.SupportsIndex):
            return self._seq[key]
        elif isinstance(key, slice):
            seq = self._seq[key]
            return type(self)(seq)
        else:
            raise TypeError(f"can't index with {type(key)}")

    def __len__(self) -> int:
        return len(self._seq)

    def __iter__(self) -> ty.Iterator[str]:
        return (chr(i) for i in self._seq)

    def __add__(self: T, other: ty.Union[str, bytes, T]) -> T:
        if isinstance(other, (str, bytes)):
            other = ensure_bytes(other)
        elif isinstance(other, type(self)):
            other = other.seq
        else:
            # prevent, e.g., concat-ing a DnaSequence to a ProteinSequence
            raise TypeError(f"can't concat {type(self)} to {type(other)}")

        return type(self)(self.seq + other)

    def __mul__(self: T, by: ty.SupportsIndex) -> T:
        return type(self)(self.seq * by)

    def __repr__(self) -> str:
        if len(self._seq) > 30:
            seq = self._seq[:30].decode("ascii") + "..."
        else:
            seq = self._seq.decode("ascii")

        return f"{self.__class__.__name__}(seq={seq!r})"

    def __str__(self) -> str:
        return self._seq.decode("ascii")

    def __bytes__(self) -> bytes:
        return self._seq

    def __eq__(self, other) -> bool:
        if isinstance(other, type(self)):
            return other._seq == self._seq
        else:
            return NotImplemented

    def __hash__(self) -> int:
        return hash((type(self), self._seq))


class ProteinSequence(BaseSequence):
    """
    A protein sequence is a sequence of IUPAC amino acid code ASCII bytes.
    """

    pass


class DnaSequence(BaseSequence):
    """
    A DNA sequence is a sequence of A, T, C, G, or IUPAC nucleotide ambiguity code ASCII bytes.
    """

    def translate(self, table: int = 1, strict: bool = False) -> ProteinSequence:
        """
        Translate a DNA sequence into a protein sequence, using the specified
        NCBI table ID.

        Raises ValueError if the table argument is invalid or any characters in
        this sequence are invalid nucleotides.

        If `strict` is true, then the input must be all `ATCG`, with no
        ambiguous nucleotides. Otherwise, a ValueError is raised.
        """

        if strict:
            seq = _translate_strict(table, self._seq)
        else:
            seq = _translate(table, self._seq)

        return ProteinSequence(seq)

    def translate_self_frames(self, table: int = 1, strict: bool = False) -> ty.List[ProteinSequence]:
        """
        Translate this DNA sequence into up to 3 protein sequences, one for each possible
        reading frame on this sense.

        May return less than 3 proteins for too-short sequences.
        For example, a sequence of length 4 only has 2 reading frames,
        and a sequence of length 2 has none.

        Can raise ValueError, see `self.translate()`

        If `strict` is true, then the input must be all `ATCG`, with no
        ambiguous nucleotides. Otherwise, a ValueError is raised.
        """

        if len(self) >= 5:
            return [
                self.translate(table, strict),
                self[1:].translate(table, strict),
                self[2:].translate(table, strict),
            ]
        elif len(self) == 4:
            return [
                self.translate(table, strict),
                self[1:].translate(table, strict),
            ]
        elif len(self) == 3:
            return [
                self.translate(table, strict),
            ]
        else:
            return []

    def translate_all_frames(self, table: int = 1, strict: bool = False) -> ty.List[ProteinSequence]:
        """
        Translate this DNA sequence into at most 6 protein sequences, one for each possible
        reading frame on this sense and the reverse complement.

        May return less than 6 proteins for too-short sequences.
        For example, a sequence of length 4 only has 2 reading frames,
        and a sequence of length 2 has none.

        Can raise ValueError, see `self.translate()`

        If `strict` is true, then the input must be all `ATCG`, with no
        ambiguous nucleotides. Otherwise, a ValueError is raised.
        """

        return [
            *self.translate_self_frames(table=table, strict=strict),
            *self.reverse_complement().translate_self_frames(table=table, strict=strict),
        ]

    def reverse_complement(self, strict: bool = False) -> "DnaSequence":
        """
        Takes the reverse complement of a DNA sequence.

        Raises ValueError if any character in this sequence is an invalid nucleotide.

        If `strict` is true, then the input must be all `ATCG`, with no
        ambiguous nucleotides. Otherwise, a ValueError is raised.
        """

        if strict:
            seq = _reverse_complement_strict(self._seq)
        else:
            seq = _reverse_complement(self._seq)
        return DnaSequence(seq)


__all__ = ["BaseSequence", "DnaSequence", "ProteinSequence"]
