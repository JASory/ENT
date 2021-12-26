//!A simple library for elementary number theory, that currently supports all primitive integer types except u128,i128,usize, isize. (Which will be supported in the next release 0.0.1)


pub mod traits;
pub mod byte;
pub mod twobytes;
pub mod fourbytes;
pub mod eightbytes;
   
    mod fjprime32;
    mod fjprime64;
    
    use crate::traits::NumberTheory;
   





