/*
   Rand generating function 
   
   as most functions are benchmarked using rand inputs it is useful to see how long it takes to generate the data separately from the function
*/

use number_theory::NumberTheory;

fn bench_rand(){

 let mut start = std::time::Instant::now();
 let mut k = 0u32;
  for _ in 0..1_000_000{
    k = u32::rng();
  }
  let mut stop = start.elapsed();
  println!("1E+6 32-bit {:?}", stop);
  
   start = std::time::Instant::now();
 let mut k = 0u64;
  for _ in 0..1_000_000{
    k = u64::rng();
  }
   stop = start.elapsed();
    println!("1E+6 64-bit {:?}", stop);
    
    
     start = std::time::Instant::now();
 let mut k = 0u128;
  for _ in 0..1_000_000{
    k = u128::rng();
  }
   stop = start.elapsed();
    println!("1E+6 128-bit {:?}", stop);
  
  
}


fn rand_prime(){
   let mut start = std::time::Instant::now();
 let mut k = 0u16;
  for _ in 0..1_000_000{
    k = u16::prime_gen(16);
  }
  let mut stop = start.elapsed();
  println!(" 1E+6 16-bit primes {:?}", stop);
  
  start = std::time::Instant::now();
   let mut k = 0u32;
  for _ in 0..1_000_000{
    k = u32::prime_gen(32);
  }
   stop = start.elapsed();
  println!(" 1E+6 32-bit primes {:?}", stop);
  
  start = std::time::Instant::now();
   let mut k = 0u64;
  for _ in 0..1_000_000{
    k = u64::prime_gen(64);
  }
   stop = start.elapsed();
  println!(" 1E+6 64-bit primes {:?}", stop);
  
  start = std::time::Instant::now();
   let mut k = 0u128;
  for _ in 0..1_000_000{
    k = u128::prime_gen(128);
  }
   stop = start.elapsed();
  println!(" 1E+6 128-bit primes {:?}", stop);
}


fn main(){
  bench_rand();
  rand_prime()
}
