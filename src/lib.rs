#![allow(non_snake_case)]
// Copyright 2021-2024 SecureDNA Stiftung (SecureDNA Foundation) <licensing@securedna.org>
// SPDX-License-Identifier: MIT OR Apache-2.0

extern crate core;

mod errors;
mod nucleotide;
pub mod trans_table; // needs to be public for bin/gen_table

pub mod canonical;

mod extendable;
pub use extendable::*;

pub mod expansions;

mod fasta;
pub use fasta::*;

mod iter;
pub use iter::*;

mod rust_api;
pub use rust_api::*;

#[cfg(feature = "python-support")]
mod python_api;

#[cfg(any(feature = "quickcheck", test))]
mod quickcheck;

#[cfg(feature = "serde")]
mod serde_utils;
