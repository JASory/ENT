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

 use std::str;
 
 fn note(x:&str,y:&str){
   use std::fs::File;
   use std::io::Write;
  let mut file = File::create(y).expect("Whoopsie, ran out of paper . . .");

   println!("My pen is ready . . . ");
   file.write_all(x.as_bytes()).expect("Whoopsie ran out of ink . . .");
}
 
 fn sqrt_precision(x: Mpz, places: u64, linelen: usize) -> String {
      let lead = x.sqrt().to_string();
      let time = std::time::Instant::now();
       let dec = Mpz::from_u64(5).pow(places*2).shl((places*2) as usize);
       let stop = time.elapsed();
      // println!("{:?}", stop);
      let tmp = x.product(dec).sqrt().to_string(); 
      let len = tmp.len();
      let trailing = tmp.as_bytes();

           let mut result : Vec<u8> = vec![];
           for i in lead.as_bytes(){
             result.push(*i)
           }
           result.push(46); // inserts decimal point  46  10
           let len = trailing.len();
           for i in 0..len{//.as_bytes()[lead.len()..].iter(){
            if i%linelen == 0 && i > 0{
              result.push(10)
            }
             result.push(trailing[i])
           }
           
           String::from_utf8(result).unwrap()
 }
 
  // 203280221primes counted in 565.992596603s
  fn main() { //1211.754987777s                ln for 1 million digits in 142.926623269s

      let time = std::time::Instant::now();
    //note(&sqrt_precision(Mpz::from_u64(5),1_000_000,40), "/home/jasory/Projects/sqrt5_2");
    let z = Mpz::new(Sign::Positive, vec![u64::MAX;2]).is_prime();
           let stop = time.elapsed();
    println!("{} {:?}",z, stop)
    
    }
