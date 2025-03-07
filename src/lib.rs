pub mod ntrait;
pub mod primitive;
pub mod data;
pub mod result;
pub mod structs;
mod arithmetic;

pub use ntrait::NumberTheory;
pub use result::NTResult;
pub use crate::arithmetic::{mpz::Mpz,sign::Sign};
