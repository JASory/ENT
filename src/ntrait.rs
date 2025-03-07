use crate::result::NTResult;
use crate::structs::{Factorization,Certificate};

pub(crate) trait MontArith{

  fn to_mont(self,ring: Self) -> Self;
  
  fn to_z(self, inv: Self, ring: Self) -> Self; 
  
  fn inv_2(self) -> Self;
  
  fn n_identity(self) -> Self;
  
  fn mont_sub(self,subtrahend: Self, ring: Self) -> Self;
  
  fn mont_add(self,addend: Self, ring: Self) -> Self;
  
  fn mont_prod(self, multiplicand: Self, inv: Self, ring: Self) -> Self;
  
  fn mont_sqr(self, inv: Self, ring: Self) -> Self;
  
  fn mont_pow(self,one: Self,p: Self, inv: Self, ring: Self) -> Self;
  
  fn odd_pow(self, p: Self, ring: Self) -> Self; 
  
  fn even_pow(self, p: Self, twofactor: Self) -> Self;
}

/// Trait for number-theory functions across all integer types
pub trait NumberTheory : Default + Clone + Sized + std::fmt::Display + std::cmp::Ord{
    /// Random number generation, generally high-quality
    fn rng() -> Self;

    ///Euclidean function, normally called euclidean division, returns Q,R such that Q*y + R = x
    /// # Panic
    /// y == 0
    fn euclidean_div(&self, other: Self) -> (Self, Self)
    where
        Self: Sized;

    /// Returns true if 1 or -1
    fn is_unit(&self) -> bool;
        
    /// Returns the representation of x in the ring Z\[n\]     
    fn residue(&self, ring: Self) -> Self;
    
    /// Returns x^-1 := x^-1*x mod n = 1 
    /// # DNE
    /// Multiplicative inverse doesn't exist
    fn mul_inverse(&self,n: Self)  -> NTResult<Self>;

    /// base^x-1 mod x == 1
    fn fermat(&self, base: Self) -> bool;
    /** Strong probable prime test.Performs Selfridge test exactly which means that if gcd(x,base) > 1 and x is a prime 
    then it will be flagged as composite

    Correct primality tests account for this, so it is up to the user to account for them as well.
    */
    fn strong_fermat(&self, base: Self) -> bool;
    /**
    Optimized for the average case. Although speed is valued over
    determinism, the probability of failure is extremely small and never deterministic beyond 2^128 (i.e repeating the test will
    almost certainly fail a previous composite that passed).

     N < 2^64 + 2^49 Provably correct deterministic test, uniquely uses a maximum of 2 strong fermat tests giving it the lowest average
     complexity publicly known. Average complexity 0.3  tests, worst-case 2.6 sprp checks (against 2^64-59)
     
     2^64 < N < 2^128 Deterministic BPSW variant, no known counterexamples but they may exist. 

     N > 2^128  Weighted to ensure 2^-64 probability of failure against a random input.Performs a minimum of 3 sprp checks. 
     Strongest counterexamples are of the form n = p(x(p-1)+1)(y(p-1)+1) where p is prime and gcd(x,y,p) = 1 and n > 2^512, passing at 
     a rate of approximately 25%. Any other equally strong counterexamples are encouraged to be reported. Further strengthening the test 
     should be by calling [sprp_check](struct.Mpz.html#method.sprp_check) afterwards **not** by calling is_prime again. Average complexity 0.18
      worst case 12 sprp checks (typical worst case is around 3.1, 12 is the absolute worst case against a very narrow interval). 


     Mersenne : Deterministic, not computed

     Fermat : Deterministic, not computed,  assumed to not exist beyond 65537 .
    */
    /* {undocumented comment} 
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
    fn prime_proof(&self) -> Certificate<Self>
    where
        Self: Sized;

    /**
         Enumerates the primes between self and sup values. If sup < self then it enumerates up to self from sup.
         Note that while evaluation is possible for any interval it becomes infeasible with very large numbers
          (>10^6000) or very large intervals (>10^12).
    */
    fn prime_list(&self, sup: Self) -> Vec<Self>
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

    /// Factorizes 
    
    /// # InfiniteSet 
    /// x == 0
    /// # DNE
    /// x == 1
    fn factor(&self) -> NTResult<Factorization<Self>>
    where
        Self: Sized+Clone;
        
    /// Returns the integer component of at least one solution to the equation x*x = n where x in Z\[i\] (aka the Gaussian integers).
    /// When   x  < 0 the result is (sqrt(x),1) otherwise it is the (sqrt(x),0)
    fn sqrt(&self) -> (Self, Self)
    where
        Self: Sized;
        
    /// Returns the integer component of one of the solutions to the equation x*x..*x where x in Z\[i\]
    fn nth_root(&self, n: Self) -> (Self, Self)
    where
        Self: Sized;
        
    /** Returns the representation of x as a base and exponent biasing towards high exponents. E.g 81 -> 3^4 
        If no higher representation of the number exists then it will return the trivial solution of x^1
    */
    fn max_exp(&self) -> (Self, Self)
    where 
        Self: Sized;
        
    /// Binary gcd, Computes the greatest common divisor of two numbers
    fn gcd(&self, other: Self) -> Self;

    /**
     Extended Euclidean GCD algorithm.
     The behavior of this function varies depending on if the type is in Z or N.

     Z : returns the GCD and Bezout coefficients (x,y) such that gcd(a,b) = ax + by

     N : returns the GCD and Bezout coefficients such that ax + by = gcd(a,b) over Z\[lcm(a,b)\]

     Note that if gcd(a,b) != a OR b this is equivalent to ax mod b  = gcd(a,b)  and by mod a = gcd(a,b).
     In the case of gcd(a,b) =1 these coefficients are the multiplicative inverses of a over Z\[b\]
      and b over Z\[a\], respectively


    */
    fn extended_gcd(&self, other: Self) -> (Self, Self, Self)
    where
        Self: Sized;
    
    /// Determines if a gcd(a,b) == 1
    fn coprime(&self, other: Self) -> bool{
          self.gcd(other).is_unit()
	}
    
    /// Computes the least common multiple checking for overflow, and zeroes
    /// # Overflow
    /// lcm(x,y) > datatype MAX
    fn lcm(&self, other: Self) -> NTResult<Self>;

   //  Exponent of the multiplicative group Z/nZ, also called the Carmichael totient of n
   //  # Infinite
   //   x == 0 As the infinite group of Z/0Z has no exponent 
    fn exponent(&self) -> NTResult<Self>
    where
        Self: Sized + Clone;
    
    /// Counts the number of coprimes from 0 to self
    fn euler_totient(&self) -> Self;

    /// Counts the Jordan's totient
    ///# Overflow
    /// Jordan_totient(x,k) > datatype MAX
    /// # CompOverflow
    /// Computation overflowed
    fn jordan_totient(&self, k: Self) -> NTResult<Self>
    where
        Self: Sized + Clone;

    ///  Higher-order Dedekind psi
    /// # Overflow
    /// dedekind_psi(x,k) > datatype MAX
    /// # CompOverflow
    /// Computational overflow
    fn dedekind_psi(&self, k: Self) -> NTResult<Self>
    where
        Self: Sized + Clone;

    /// Returns x*y mod n
    /// # Failure
    /// if x * y > datatype MAX  AND n == 0 
    fn product_residue(&self, other: Self, n: Self) -> Self;
    
    /// Returns x*y mod n
    /// # Overflow
    /// if x * y > datatype MAX  AND n == 0 
    fn checked_product_residue(&self, other: Self, n: Self) -> NTResult<Self>
    where
        Self: Sized + Clone;
     
    /// Returns x*x mod n, similar to product_residue except more optimized squaring.
    /// # Failure
    /// if x * x > datatype MAX  AND n == 0
    fn quadratic_residue(&self, n: Self) -> Self;

    /// Returns x*x mod n, similar to product_residue except more optimized squaring.
    /// # Overflow
    /// if x * x > datatype MAX  AND n == 0
    fn checked_quadratic_residue(&self, n: Self) -> NTResult<Self>
    where
        Self: Sized + Clone;
        
    /** Returns x^y mod n, generally called mod_pow in other libraries. If y < 0 returns -x^|y| mod n,
    aka the exponentiation of the multiplicative inverse, analogous to the behavior in the Reals.
   */ 
    /// # Failure
    /// If y < 0 and gcd(x,n) > 1, or n = 0 and x^y > datatype MAX
    
    fn exp_residue(&self, pow: Self, n: Self) -> Self;

    /// Exponential residue x^y mod n
    /// # DNE 
    ///  y < 0 AND gcd(x,n) > 1
    /// # Overflow
    /// n == 0 AND x^y > datatype MAX 
    fn checked_exp_residue(&self, pow: Self, n: Self) -> NTResult<Self>
    where
        Self: Sized + Clone;

    // Returns an x such that x*x = a over Z/nZ    
    // # None
    // If x does not exist
 //   fn sqrt_residue(&self, n: Self) -> Option<Self>
 //   where 
 //       Self: Sized;
        
    /// Determines of a number is k-free, i.e square-free, cube-free etc. . .
    fn k_free(&self, k: Self) -> bool;
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
    /// # Infinite
    /// x == 0 
    fn smooth(&self) -> NTResult<Self>
    where
        Self: Sized + Clone;


    /// Checks if the smoothness bound is at least b
    fn is_smooth(&self, b: Self) -> bool;

    /// Legendre symbol of a,p.
    /// # Failure
    /// P is not an odd prime
    fn legendre(&self, p: Self) -> i8;
    
    /// Legendre symbol of a,p. 
    /// # Undefined 
    /// P is not an odd prime
    fn checked_legendre(&self, p: Self) -> NTResult<i8>;

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
    fn jacobi(&self, p: Self) -> i8;

    /// Jacobi symbol of x,p.
    /// # Undefined
    /// P is not an odd, positive integer
    fn checked_jacobi(&self, p: Self) -> NTResult<i8>;
    
    /// Kronecker symbol
    fn kronecker(&self, k: Self) -> i8;
    
    /// Multiplicative order of a mod N
    /// # DNE
    /// Multiplicative order does not exist
    fn ord(&self,n: Self) -> NTResult<Self>;
}

// ENT maps to 32-bit as modular arithmetic requires promoting to 64-bit

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

reducer!(u64;u32, i64;i32, u128;u64, i128;i64);
/*
impl Reduction for Mpz{
  fn reducible(&self) -> bool{
    if self.len() < 3{
      return true
    }
    false
  }
}
*/
/*
pub(crate) trait NumberTheory{
    fn fermat
    fn strong_fermat
    fn gcd
    fn jacobi
    fn legendre
    fn kronecker
    fn prime_proof
    fn is_prime
    fn sqrt_residue
    fn exp_residue
    fn extended_gcd
    fn nth_prime
    fn pi
    fn prime_list
    fn prime_gen
    fn factor
    fn lcm
    fn exponent
    fn euler
    fn jordan
    fn dedekind
    fn radical
    fn k_free
    fn mangoldt
    fn mobius
    fn power_residue
    fn ind
    fn ord
    fn partition
    fn primitive_root
}
*/


