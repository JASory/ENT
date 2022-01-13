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
    
    use crate::traits::NumberTheory;
    use crate::arithmetic::mpz::Mpz;
    use crate::arithmetic::sign::Sign;
   // use crate::arithmetic::inlince
  
 //247330401473104534033686901979090442377579073428791223965108330490992066367285383405085116763961027192036158782418658572232976002286069292939485950791477815115016413693758078
 






