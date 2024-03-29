/*
  Traits

*/
use crate::NTResult; 
use crate::Mpz;

/// Trait for number-theory functions across all integer types
pub trait NumberTheory : Default + Clone + Sized + std::fmt::Display{
    ///Random number generation, generally high-quality
    fn rng() -> Self;

    ///Euclidean function, normally called euclidean division, returns Q,R such that Q*y + R = x
    /// # Panic
    /// y == 0
    fn euclidean_div(&self, other: &Self) -> (Self, Self)
    where
        Self: Sized;
        
    /// Returns the representation of x in the ring Z\[n\]     
    fn residue(&self, ring: &Self) -> Self; 
       
    /** Strong probable prime test.Performs Artjuhov-Selfridge test exactly which means that if gcd(x,base) > 1 and x is a prime 
    then it will be flagged as composite

    Correct primality tests account for this, so it is up to the user to account for them as well.
    */
    fn is_sprp(&self, base: &Self) -> bool;
    /**
    Optimized for the average case. Although speed is valued over
    determinism, the probability of failure is extremely small and never deterministic (i.e repeating the test will
    almost certainly fail a previous composite that passed).

     N < 2^64 + 2^49 Provably correct deterministic test, uniquely uses a maximum of 2 strong fermat tests giving it the lowest average
     complexity publicly known. Average complexity 0.3  tests, worst-case 2.6 sprp checks (against 2^64-59)


     N > 2^64 + 2^49  Weighted to ensure 2^-64 probability of failure against a random input.Performs a minimum of 3 sprp checks. 
     Strongest counterexamples are of the form n = p(x(p-1)+1)(y(p-1)+1) where p is prime and gcd(x,y,p) = 1 and n > 2^512, passing at 
     a rate of approximately 25%. Any other equally strong counterexamples are encouraged to be reported. Further strengthening the test 
     should be by calling [sprp_check](struct.Mpz.html#method.sprp_check) afterwards **not** by calling is_prime again. Average complexity 0.18
      worst case 12 sprp checks (typical worst case is around 3.1, 12 is the absolute worst case against a very narrow interval). 


     Mersenne : Deterministic, not computed

     Fermat : Deterministic, not computed,  assumed to not exist beyond 65537 .
    */
    /* {undocumented comment} is_prime is further notable for having the following properties : 
        highest single-shot bound (2^35 or eight times higher than 2^32 which is the standard in all other hashtable tests)
        highest bound for any test with less than 12 checks (except possibly the Baille-PSW test, although it is not a MR test)
        the first correct 2-shot test (see /data/hashtable for more information on assessment and comparison with other attempts)
        Low memory relative to it's performance. 
        One might be curious as to why 3 sprp tests are selected,  this is due to a  special configuration of checks that makes them 
        equivalent to around 5 standard tests, so in reality not only does is_prime meet the advertised requirements it far exceeds it, being
        fully deterministic against nearly half of all pseudoprimes. Additionally Lucas sequence and quadratic Frobenius tests are not used
        because the current configuration is trivially parallelizable to the cost of only a single sprp test. (versus two for a lucas and
         3 for a Frobenius). As stated in the documented comment is_prime is configured to be optimal in the average case, the Carmichael 
         numbers that do pass with a relatively high probability are statistically so rare that they are inconsequential.
    */    
    fn is_prime(&self) -> bool;

    /**
       Checks primality for the integer and returns the evaluation in addition to a vector of the witness
       followed by the factors of n-1. For instance  269u64.prime_proof may return (true, [238, 2,67]).
        Verification requires checking that 238.exp_residue(269-1, 269) == 1
        and 238.exp_residue((269-1)/2, 269) != 1 and 238.exp_residue((269-1)/67, 269) != 1. See prime_proof in examples to see how.
        Unable to construct a proof for Carmichael numbers. 
    */
    fn prime_proof(&self) -> (bool, Vec<Self>)
    where
        Self: Sized;

    /**
         Enumerates the primes between self and sup values. If sup < self then it enumerates up to self from sup.
         Note that while evaluation is possible for any interval it becomes infeasible with very large numbers
          (>10^6000) or very large intervals (>10^12).
    */
    fn prime_list(&self, sup: &Self) -> Vec<Self>
    where
        Self: Sized;

    /**
    N-th prime, exact evaluation for less than 2^64, approximation beyond that.  Not currently feasible beyond 10^12
     */
     /// # DNE
     /// n == 0  (there is no zeroth prime)
     /// # Overflow
     /// Pn > datatype MAX

    fn nth_prime(&self) -> NTResult<Self>
    where
        Self: Sized + Clone;
        
        
    /// Generates an odd positive prime in the interval 2^(k-1);2^k 
    /// # DNE
    /// x == 0
    /// # Overflow
    /// x > datatype length in bits 
    fn prime_gen(k: u32) -> NTResult<Self>
    where 
        Self : Sized + Clone;
 
 
    /// Prime-counting function, exact evaluation for less than 2^64, approximation beyond that. Not currently feasible beyond 10^12
    fn pi(&self) -> Self;

    /// Factorizes into a vector of the form prime factor, power, prime factor, power . . . i.e 2,6,5,2 for 2^6 * 5^2  = 1600
    
    /// # Failure 
    /// x == 0 OR  x == 1
    fn factor(&self) -> Vec<Self>
    where
        Self: Sized;
    
    /// Factorizes into a vector of the form prime factor, power, prime factor, power . . . i.e 2,6,5,2 for 2^6 * 5^2  = 1600
    
    /// # InfiniteSet 
    /// x == 0
    /// # DNE
    /// x == 1
    fn checked_factor(&self) -> NTResult<Vec<Self>>
    where
        Self: Sized + Clone;    
    /// Returns the integer component of at least one solution to the equation x*x = n where x in Z\[i\] (aka the Gaussian integers).
    /// When   x  < 0 the result is (sqrt(x),1) otherwise it is the (sqrt(x),0)
    fn sqrt(&self) -> (Self, Self)
    where
        Self: Sized;
        
    /// Returns the integer component of one of the solutions to the equation x*x..*x where x in Z\[i\]
    fn nth_root(&self, n: &Self) -> (Self, Self)
    where
        Self: Sized;
        
    /** Returns the representation of x as a base and exponent biasing towards high exponents. E.g 81 -> 3^4 
        If no higher representation of the number exists then it will return the trivial solution of x^1
    */
    fn max_exp(&self) -> (Self, Self)
    where 
        Self: Sized;
    /// Binary gcd, Computes the greatest common divisor of two numbers
    fn gcd(&self, other: &Self) -> Self;

    /**
     Extended Euclidean GCD algorithm.
     The behavior of this function varies depending on if the type is in Z or N.

     Z : returns the GCD and Bezout coefficients (x,y) such that gcd(a,b) = ax + by

     N : returns the GCD and Bezout coefficients such that ax + by = gcd(a,b) over Z\[lcm(a,b)\]

     Note that if gcd(a,b) != a OR b this is equivalent to ax mod b  = gcd(a,b)  and by mod a = gcd(a,b).
     In the case of gcd(a,b) =1 these coefficients are the multiplicative inverses of a over Z\[b\]
      and b over Z\[a\], respectively


    */
    fn extended_gcd(&self, other: &Self) -> (Self, Self, Self)
    where
        Self: Sized;

    ///  Computes the least common multiple of x,y
   ///  # Failure
   ///   If lcm(x,y) >  datatype MAX
   ///
    fn lcm(&self, other: &Self) -> Self;

    /// Computes the least common multiple checking for overflow, and zeroes
    /// # Overflow
    /// lcm(x,y) > datatype MAX
    fn checked_lcm(&self, other: &Self) -> NTResult<Self>
    where
        Self: Sized + Clone;
        
   ///  Carmichael totient function, also the exponent of the multiplicative group Z/nZ
   ///  # Infinite
   ///   x == 0 As the infinite group of Z/0Z has no exponent 
    fn carmichael_totient(&self) -> NTResult<Self>
    where
        Self: Sized + Clone;
    
    /// Counts the number of coprimes from 0 to self
    fn euler_totient(&self) -> Self;

    /// Counts the Jordan's totient
    ///# Overflow
    /// Jordan_totient(x,k) > datatype MAX
    /// # CompOverflow
    /// Computation overflowed
    fn jordan_totient(&self, k: &Self) -> NTResult<Self>
    where
        Self: Sized + Clone;

    ///  Higher-order Dedekind psi
    /// # Overflow
    /// dedekind_psi(x,k) > datatype MAX
    /// # CompOverflow
    /// Computational overflow
    fn dedekind_psi(&self, k: &Self) -> NTResult<Self>
    where
        Self: Sized + Clone;

    /// Returns x*y mod n
    /// # Failure
    /// if x * y > datatype MAX  AND n == 0 
    fn product_residue(&self, other: &Self, n: &Self) -> Self;
    
    /// Returns x*y mod n
    /// # Overflow
    /// if x * y > datatype MAX  AND n == 0 
    fn checked_product_residue(&self, other: &Self, n: &Self) -> NTResult<Self>
    where
        Self: Sized + Clone;
     
    /// Returns x*x mod n, similar to product_residue except more optimized squaring.
    /// # Failure
    /// if x * x > datatype MAX  AND n == 0
    fn quadratic_residue(&self, n: &Self) -> Self;

    /// Returns x*x mod n, similar to product_residue except more optimized squaring.
    /// # Overflow
    /// if x * x > datatype MAX  AND n == 0
    fn checked_quadratic_residue(&self, n: &Self) -> NTResult<Self>
    where
        Self: Sized + Clone;
        
    /** Returns x^y mod n, generally called mod_pow in other libraries. If y < 0 returns -x^|y| mod n,
    aka the exponentiation of the multiplicative inverse, analogous to the behavior in the Reals.
   */ 
    /// # Failure
    /// If y < 0 and gcd(x,n) > 1, or n = 0 and x^y > datatype MAX
    
    fn exp_residue(&self, pow: &Self, n: &Self) -> Self;

    /// Exponential residue x^y mod n
    /// # DNE 
    ///  y < 0 AND gcd(x,n) > 1
    /// # Overflow
    /// n == 0 AND x^y > datatype MAX 
    fn checked_exp_residue(&self, pow: &Self, n: &Self) -> NTResult<Self>
    where
        Self: Sized + Clone;
        /*
    /// Returns an x such that x*x = a over Z/nZ    
    /// # None
    /// If x does not exist
    fn sqrt_residue(&self, n: &Self) -> Option<Self>
    where 
        Self: Sized;
        */
    /// Determines of a number is k-free, i.e square-free, cube-free etc. . .
    fn k_free(&self, k: &Self) -> bool;
/*
    fn partition(&self) -> Option<Self>
    where
        Self: Sized;
  */
    
    /// Returns the product of each prime such that p | n AND p > 0
    /// # Infinite
    ///  x == 0
    fn radical(&self) -> NTResult<Self>
    where
        Self: Sized + Clone;

    /// Returns the smoothness bound of n, this is the largest prime factor of n. 
    fn smooth(&self) -> NTResult<Self>
    where
        Self: Sized + Clone;
    
    //fn checked_smooth(&self)

    /// Checks if the smoothness bound is at least b
    fn is_smooth(&self, b: &Self) -> bool;

    /// Legendre symbol of a,p.
    /// # Failure
    /// P is not an odd prime
    fn legendre(&self, p: &Self) -> i8;
    
    /// Legendre symbol of a,p. 
    /// # Undefined 
    /// P is not an odd prime
    fn checked_legendre(&self, p: &Self) -> NTResult<i8>;

    /// Liouville function
    fn liouville(&self) -> i8;
    
    /// Lagarias derivative
    /// # Overflow
    /// D(x) > datatype MAX
    fn derivative(&self) -> NTResult<Self>
    where 
        Self : Sized + Clone;

    ///Von Mangoldt function, returns the natural logarithm of self if it is a prime-power
    fn mangoldt(&self) -> f64;

    ///  Mobius function
    fn mobius(&self) -> i8;

    /// Jacobi symbol of x,p. 
    /// # Failure
    /// P is not an odd, positive integer
    fn jacobi(&self, p: &Self) -> i8;

    /// Jacobi symbol of x,p.
    /// # Undefined
    /// P is not an odd, positive integer
    fn checked_jacobi(&self, p: &Self) -> NTResult<i8>;
    
    /// kronecker symbol
    fn kronecker(&self, k: &Self) -> i8;
    /*
    /// Pi approximation
    fn pi(&self) -> f64;

    let x = self.ln();
     //(self/x )*(1.0 + 1.0/x  + 2.0/(x.ln()*x.ln()))

    //Primitive root


    Future functions


kronecker symbol
power residue
partition
multiplicative order
primitive root
primitive root of unity
index
sqrt residue
divisor list (list of all proper divisors)


index, kronecker, partition, mul_order, sqrt
    */
}

//const trait hello{}

pub(crate) trait Reduction{
  fn reducible(&self) -> bool;
}

macro_rules! reducer(
  ($($t:ty; $s:ty),* $(,)*) => {$(
  
  impl Reduction for $t{
  #[inline]
  fn reducible(&self) -> bool{
    if *self < <$s>::MAX as $t{
      return true
    }
    return false
  } 
  }
    )*}
);

reducer!(u16;u8, i16;i8, u32;u16, i32;i16, u64;u32, i64;i32, u128;u64, i128;i64 );

impl Reduction for Mpz{
  fn reducible(&self) -> bool{
    if self.len() < 3{
      return true
    }
    false
  }
}

/*
NTT research

  https://stackoverflow.com/questions/18465326/fast-bignum-square-computation

  https://stackoverflow.com/questions/18577076/modular-arithmetics-and-ntt-finite-field-dft-optimizations
  https://stackoverflow.com/questions/10260927/translation-from-complex-fft-to-finite-field-fft/18547575#18547575

  https://dl.acm.org/doi/10.1145/1185448.1185568

*/
