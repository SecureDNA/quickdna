# Copyright 2021-2024 SecureDNA Stiftung (SecureDNA Foundation) <licensing@securedna.org>
# SPDX-License-Identifier: MIT OR Apache-2.0

import sys

from Bio.Seq import Seq
import timeit

sys.path.insert(0, ".")
from quickdna import DnaSequence

with open("benchmarks/covid.txt", "rb") as f:
    covid_genome = f.read().replace(b"\n", b"")

small_genome = covid_genome[:300]


def translate_biopython(dna):
    _translated = Seq(dna).translate(table=1)


def translate_quickdna(dna):
    _translated = DnaSequence(dna).translate(table=1)


def reverse_complement_biopython(dna):
    _rc = Seq(dna).reverse_complement()


def reverse_complement_quickdna(dna):
    _rc = DnaSequence(dna).reverse_complement()


ITERS = 1_000


def report_bench(callback_name: str, argument_name: str, baseline=None) -> float:
    code = (
        f"from __main__ import {callback_name}, {argument_name};\n"
        f"{callback_name}({argument_name})"
    )

    time = timeit.timeit(code, number=ITERS) / ITERS * 1_000
    if baseline:
        comparison = f" ({time / baseline * 100:.2f}%)"
    else:
        comparison = ""

    print(f"{callback_name}({argument_name}) = {time:.5f}ms / iter{comparison}")
    return time


def report_bench_for(callback_prefix: str, argument_name: str, do_compare: bool):
    baseline = report_bench(f"{callback_prefix}_quickdna", argument_name)
    if do_compare:
        report_bench(f"{callback_prefix}_biopython", argument_name, baseline)


if __name__ == "__main__":
    do_compare = "--no-compare" not in sys.argv

    report_bench_for("translate", "small_genome", do_compare)
    report_bench_for("translate", "covid_genome", do_compare)
    report_bench_for("reverse_complement", "small_genome", do_compare)
    report_bench_for("reverse_complement", "covid_genome", do_compare)
