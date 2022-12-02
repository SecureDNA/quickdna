#![allow(non_snake_case)]

extern crate core;

mod errors;
mod nucleotide;
pub mod trans_table; // needs to be public for bin/gen_table

mod fasta;
pub use fasta::*;

mod iter;
pub use iter::*;

mod rust_api;
pub use rust_api::*;

#[cfg(feature = "python-support")]
mod python_api;

#[cfg(feature = "python-support")]
pub use python_api::*;
