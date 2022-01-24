//!A simple library for elementary number theory and arbitrary precision arithmetic. 


pub mod traits;
pub mod byte;
pub mod twobytes;
pub mod fourbytes;
pub mod eightbytes;
pub mod sixteenbytes;   
    mod fjprime32;
    mod fjprime64;
    mod primes;
    
    // computational speed 193.483187015s for primality checks in the interval [0;10^9]
 pub mod arithmetic;
    
 pub   use crate::traits::NumberTheory;
 pub   use crate::arithmetic::mpz::Mpz;
 pub   use crate::arithmetic::sign::Sign;

  // 203280221primes counted in 565.992596603s
  

