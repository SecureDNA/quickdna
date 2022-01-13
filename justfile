# requires `just`: cargo install just # https://github.com/casey/just

default:
    just --list

alias t := test
alias c := check

# Build a wheel with maturin
build:
    poetry run maturin build --release

# Builds a wheel with maturin and installs it to the current virtualenv
develop:
    poetry run maturin develop --release

# Builds the wheel, installs it, and runs pytest
test: develop
    HYPOTHESIS_PROFILE=dev pytest

# Builds the wheel, installs it, and runs the benchmarks against biopython
bench: develop
    python benchmarks/bench.py

# Builds the wheel, installs it, and runs cargo-flamegraph to generate flamegraph.svg
profile: develop
    flamegraph python benchmarks/bench.py --no-compare

# Runs fmt, clippy, and black
check: fmt clippy black

# Runs cargo fmt to format Rust code
fmt:
    cargo fmt --all

# Runs clippy to lint Rust code
clippy:
    cargo clippy --all

# Runs black to format Python code
black:
    poetry run black .
