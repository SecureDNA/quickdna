import typing as ty

from .quickdna import _translate, _reverse_complement  # type: ignore

T = ty.TypeVar("T", bound="BaseSequence")

class BaseSequence:
    """Base class for DNA and Protein sequences."""

    def __init__(self, seq: ty.Union[str, bytes]) -> None:
        if isinstance(seq, str):
            seq = seq.encode("ascii", errors="strict")
        self._seq = seq

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
            raise ValueError("unknown __getitem__ key", key)

    def __repr__(self) -> str:
        if len(self._seq) > 30:
            seq = self._seq[:30].decode("ascii") + "..."
        else:
            seq = self._seq.decode("ascii")

        return f"{self.__class__.__name__}(seq={seq!r})"

    def __str__(self) -> str:
        return self._seq.decode("ascii")

    def __eq__(self, other) -> bool:
        if isinstance(other, type(self)):
            return other._seq == self._seq
        else:
            return NotImplemented

    def __hash__(self) -> int:
        return hash((type(self), self._seq))


class ProteinSequence(BaseSequence):
    pass


class DnaSequence(BaseSequence):
    def translate(self, table: int = 1) -> ProteinSequence:
        """
        Translate a DNA sequence into a protein sequence, using the specified
        NCBI table ID.
        """

        seq = _translate(table, self._seq)
        return ProteinSequence(seq)

    def translate_all_frames(
        self, table: int = 1
    ) -> ty.Tuple[ProteinSequence, ProteinSequence, ProteinSequence]:
        """
        Translate this DNA sequence into 3 protein sequences, one for each possible
        reading frame.
        """

        return (
            self.translate(table),
            self[1:].translate(table),
            self[2:].translate(table),
        )

    def reverse_complement(self) -> "DnaSequence":
        """Takes the reverse complement of a DNA sequence"""

        seq = _reverse_complement(self._seq)
        return DnaSequence(seq)


__all__ = ["BaseSequence", "DnaSequence", "ProteinSequence"]
