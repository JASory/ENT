/*
  Traits

*/

/// Trait for number-theory functions across all integer types
pub trait NumberTheory {
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

     N < 2^64 + 2^47 Provably Deterministic, uniquely uses a maximum of 2 strong fermat tests giving it the lowest average complexity 
    publicly known. Average complexity 0.3  tests, worst-case 2.6 sprp checks (against 2^64-59)


     N > 2^64 + 2^47  Weighted to ensure 2^-64 probability of failure against a random input.Performs a minimum of 3 sprp checks. 
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
        Unable to construct a proof for the Carmichael numbers.
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
     /// # None
     /// Pn > datatype MAX

    fn nth_prime(&self) -> Option<Self>
    where
        Self: Sized;
    
    /*
    fn checked_prime_gen(k: u32) -> Option<Self>
    where 
        Self: Sized;
    */    
        
    /// Generates an odd positive prime in the interval 2^(k-1);2^k . 
    fn prime_gen(k: u32) -> Self;
 
    /// Prime-counting function, exact evaluation for less than 2^64, approximation beyond that. Not currently feasible beyond 10^12
    fn pi(&self) -> Self;

    /// Factorizes into a vector of the form factor, power, factor, power . . . i.e 2,6,5,2 for 2^6 * 5^2  = 1600
    fn factor(&self) -> Vec<Self>
    where
        Self: Sized;
        
    /// Returns the integer component of at least one solution to the equation x*x = n where x in Z\[i\] (aka the Gaussian integers).
    /// When   x  < 0 the result is (sqrt(x),1) otherwise it is the (sqrt(x),0)
    fn sqrt(&self) -> (Self, Self)
    where
        Self: Sized;
        
    /// Returns the integer component of one of the solutions to the equation x*x..*x where x in Z\[i\]
    fn nth_root(&self, n: &Self) -> (Self, Self)
    where
        Self: Sized;

    /// Binary gcd, Computes the greatest common divisor of two numbers
    fn euclid_gcd(&self, other: &Self) -> Self;

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
    /// # None
    /// lcm(x,y) > datatype MAX
    fn checked_lcm(&self, other: &Self) -> Option<Self>
    where
        Self: Sized;
        
    // Carmichael function 
   // fn carmichael(&self) -> Self;
    
    /// Counts the number of coprimes from 0 to self
    fn euler_totient(&self) -> Self;

    /// Counts the Jordan's totient
    ///# None 
    /// Jordan_totient(x,k) > datatype MAX
    fn jordan_totient(&self, k: &Self) -> Option<Self>
    where
        Self: Sized;

    ///  Higher-order Dedekind psi
    /// # None
    /// dedekind_psi(x,k) > datatype MAX
    fn dedekind_psi(&self, k: &Self) -> Option<Self>
    where
        Self: Sized;

    /// Returns x*y mod n

    fn product_residue(&self, other: &Self, n: &Self) -> Self;

    /// Returns x*x mod n, similar to product_residue except more optimized squaring.
    fn quadratic_residue(&self, n: &Self) -> Self;

    /** Returns x^y mod n, generally called mod_pow in other libraries. If y < 0 returns -x^|y| mod n,
    aka the exponentiation of the multiplicative inverse, analogous to the behavior in the Reals.
   */ 
    /// # Failure
    /// If y < 0 and gcd(x,n) > 1, or n = 0 and x^y > datatype MAX
    
    fn exp_residue(&self, pow: &Self, n: &Self) -> Self;

    /// Exponential residue checked for overflow (in the case of Z/0Z) and lack of an inverse x^-1 in the case of y < 0
    fn checked_exp_residue(&self, pow: &Self, n: &Self) -> Option<Self>
    where
        Self: Sized;

    /// Determines of a number is k-free, i.e square-free, cube-free etc. . .
    fn k_free(&self, k: &Self) -> bool;
/*
    fn partition(&self) -> Option<Self>
    where
        Self: Sized;
        */
    /// Returns the product of each prime such that p | n
    fn radical(&self) -> Self;

    /// Returns the smoothness bound of n, this is the largest prime factor of n. 
    fn smooth(&self) -> Self;

    /// Checks if the smoothness bound is at least b
    fn is_smooth(&self, b: &Self) -> bool;

    /// Legendre symbol of a,p.
    /// # Failure
    /// P is not an odd prime
    fn legendre(&self, p: &Self) -> i8;
    
    /// Legendre symbol of a,p. 
    /// # None 
    /// P is not an odd prime
    fn checked_legendre(&self, p: &Self) -> Option<i8>;

    /// Liouville function
    fn liouville(&self) -> i8;
    
    /// Lagarias derivative
    fn derivative(&self) -> Option<Self>
    where 
        Self : Sized;

    ///Von Mangoldt function, returns the natural logarithm of self if it is a prime-power
    fn mangoldt(&self) -> f64;

    ///  Mobius function
    fn mobius(&self) -> i8;

    /// Jacobi symbol of x,p. 
    /// # Failure
    /// P is not an odd, positive integer
    fn jacobi(&self, p: &Self) -> i8;

    /// Jacobi symbol of x,p.
    /// # None
    /// P is not an odd, positive integer
    fn checked_jacobi(&self, p: &Self) -> Option<i8>;
    /*

    /// Pi approximation
    fn pi(&self) -> f64;

    let x = self.ln();
     //(self/x )*(1.0 + 1.0/x  + 2.0/(x.ln()*x.ln()))

    //Primitive root


    Future functions


kronecker symbol
power residue
carmichael
partition
multiplicative order
primitive root
primitive root of unity
index
sqrt residue


    */
}

/*
NTT research

  https://stackoverflow.com/questions/18465326/fast-bignum-square-computation

  https://stackoverflow.com/questions/18577076/modular-arithmetics-and-ntt-finite-field-dft-optimizations
  https://stackoverflow.com/questions/10260927/translation-from-complex-fft-to-finite-field-fft/18547575#18547575

  https://dl.acm.org/doi/10.1145/1185448.1185568

*/
