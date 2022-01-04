//!A simple library for elementary number theory and 


pub mod traits;
pub mod byte;
pub mod twobytes;
pub mod fourbytes;
pub mod eightbytes;
pub mod sixteenbytes;   
    mod fjprime32;
    mod fjprime64;
    mod primes;
    
    // computational speed 193.483187015s for 10^9
 pub mod arithmetic;
    
    use crate::traits::NumberTheory;
    use crate::arithmetic::mpz::Mpz;
    use crate::arithmetic::sign::Sign;
    
    fn main(){
   
  
  
  println!("{:?}",Mpz::sirp(1,10000,2,0).to_string());//22207  22475  22478  22489
    }
   





