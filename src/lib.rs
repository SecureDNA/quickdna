#![feature(array_chunks)]
#![feature(core_intrinsics)]
#![allow(non_snake_case)]

mod errors;
mod nucleotide;
pub mod trans_table; // needs to be public for bin/gen_table

mod rust_api;
pub use rust_api::*;

mod python_api;
pub use python_api::*;
