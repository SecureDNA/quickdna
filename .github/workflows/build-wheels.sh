#!/bin/bash
set -euxo pipefail nullglob

for PY in cp38-cp38 cp39-cp39 cp310-cp310; do
    PYBIN="/opt/python/${PY}/bin"
    "${PYBIN}/pip" install maturin
    "${PYBIN}/maturin" build -i "${PYBIN}/python" --release
done

for wheel in target/wheels/*.whl; do
    auditwheel repair "${wheel}"
done
