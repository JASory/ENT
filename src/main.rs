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
 


fn hardware_rng()-> u64{
    let mut x: u64 = 0; 
    let k = unsafe { core::arch::x86_64::_rdrand64_step(&mut x) } ;
   x
 }
  
  
  fn main(){
  // 203280221primes counted in 900.660736901s // 203280221 203152305.728129     999897657631
  // 203280221primes counted in 549.946332932s
  //4506732primes counted in 38.294082118s

  // 2^86243-1 in  2276.469896756s
  let t = Mpz::from_u64(10);
   let mut count = 0u64;
   let mut k = Mpz::new(Sign::Positive, vec![u64::MAX;2]);
   k.mut_addition(t);
   //println!("{}",k.to_string());
   
   
let time = std::time::Instant::now();
for i in 0..300_000{
 if k.probable_prime(){
  count+=1;
  //println!("13^ 8258 + {} is prime",i+10);
 }
 //println!("Evaluated 13^ 8258 + {}",i);
 k.successor();
}
let stop = time.elapsed();
println!("{:?}primes counted in {:?}",count, stop)
  
  }





