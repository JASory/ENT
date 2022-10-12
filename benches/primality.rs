use number_theory::NumberTheory;
use number_theory::Mpz;

/*
Benchmark competition and show errors
  Rivals
  num-modular  legendre error  (Cmpute's libraries are completely (but silently) broken)
  num-prime    is_prime_error   32k+ counterexamples
  num-primes   prime error ?  (Not a typo, num-primes is an older attempt  by  a different author at a simple primality checking crate, 
                                it too was broken although it's current status has not been evaluated)
  num-integer  Compare gcd 
  num-bigint    
  malachite                 The only real comparable library
  primal     is_prime error (since fixed due to the author)
  primetools
  red_primality
  glass-pumpkin      Baille-psw (lucas test) is broken and the miller rabin over small intervals is far too weak (atleast in old versions) 

*/
// 64-bit known composites that deterministically pass some of the libraries, this is not an exhaustive list (some of them have thousands to millions of known counterexamples )
//const COUNTEREXAMPLES : [u64;]

/* An oft looked over property is that compositeness tests like miller-rabin and lucas sequence tests, can infact flag primes as composite
 especially if they are not correctly implemented. this is harder to detect than passing composites however and is primarily found by code auditing
  rather than running against a list of primes*. glass-pumpkin (1.2) BPSW for instance flags primes due to an incorrect lucas test, as will any miller-rabin check that fails to guarantee that a prime is only checked by a coprime base. 
  *Of course a fast provably correct library like number-theory can also help detect these errors, and even better is that number-theory can also produce a proof (via prime_proof function) even if you don't trust the is_prime function. 
   
*/
pub fn bench_prime() {
    let mut count = 0i64;
    println!("Evaluating 1 million integers in the worst interval for each datatype");
    let mut start = std::time::Instant::now();
    for i in 4293967295..u32::MAX{
        if i.is_prime() {
            count += 1
        }
    }
    let mut stop = start.elapsed();
    println!(
        "{} primes counted in {:?} with an error of {} : u32",
        count,
        stop,
        count - 44872
    );
    count = 0;
    start = std::time::Instant::now();
    for i in 18446744073708551615..u64::MAX {
        if i.is_prime() {
            count += 1
        }
    }
    stop = start.elapsed();
    println!(
        "{} primes counted in {:?} with an error of {} : u64",
        count,
        stop,
        count - 22475
    );

    count = 0;
    start = std::time::Instant::now();
    for i in (u128::MAX - 1_000_000)..u128::MAX {
        if i.is_prime() {
            count += 1
        }
    }
    stop = start.elapsed();
    println!(
        "{} primes counted in {:?} with an error of {} : u128",
        count,
        stop,
        count - 11281
    );
    
    let mut p = Mpz::one().shl(256).ref_subtraction(&Mpz::from_u64(1_000_000));
    count = 0;
        start = std::time::Instant::now();
    for i in  0..1_000_000{
      if p.is_prime(){
        count+=1;
      }
      p.successor()
    }
    stop = start.elapsed();
    println!("{} primes counted in  {:?} with an error of {} : 256-bit Mpz",count,stop, count - 5539);
    
    let mut p = Mpz::from_u64(10).pow(999);
    
    count = 0;
            start = std::time::Instant::now();
            for i in  0..10_000{
             if p.is_prime(){
               count+=1;
             }
              p.successor()
           }
           let stop = start.elapsed();
           println!("{} titanic primes found in {:?}",count, stop)
}

fn main() {
    bench_prime()
}
