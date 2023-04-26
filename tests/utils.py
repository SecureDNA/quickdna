# Copyright 2021-2023 SecureDNA Stiftung (SecureDNA Foundation) <license@securedna.org>
# SPDX-License-Identifier: MIT OR Apache-2.0

import typing as ty


def valid_tables() -> ty.Generator[int, None, None]:
    for table in range(1, 34):
        # skip tables that are either unassigned (17-20) or are aliases for other
        # tables that biopython doesn't understand (7 & 8)
        if table not in (7, 8, 17, 18, 19, 20):
            yield table
