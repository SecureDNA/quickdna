# Copyright 2021-2023 SecureDNA Stiftung (SecureDNA Foundation) <licensing@securedna.org>
# SPDX-License-Identifier: MIT OR Apache-2.0

set -euxo pipefail nullglob

for PY in cp310-cp310 cp311-cp311 cp312-cp312; do
    PYBIN="/opt/python/${PY}/bin"
    "${PYBIN}/pip" install maturin
    "${PYBIN}/maturin" build -i "${PYBIN}/python" --release
done

for wheel in target/wheels/*.whl; do
    auditwheel repair "${wheel}"
done
