/*
  Traits

*/

pub trait NumberTheory{
   
   ///Random number generation, generally high-quality
   fn rng() -> Self;
   /// Euclidean function, normally called euclidean division
   fn euclidean_div(&self, other: &Self) -> (Self,Self) where Self: Sized;
   ///Returns optimally fast primality check for each interval. 
   /// Deterministic for n < 2^64. For n > 2^64 then the probability of failure is > 2^-64. Unlike most libraries Number-Theory weights the number of Strong Fermat checks
   /// to guarantee a constant error rate, resulting in being slightly slower, but more accurate, for small numbers but faster for larger numbers. Values greater than 2^512
   /// receive only a Strong Fermat base-2 check and a random base-check as probabilistically only a single random Strong Fermat test is needed to have an error rate of 2^-64

   fn is_prime(&self) -> bool;
   
   /// Factorizes into a vector of the form factor, power, factor, power . . . i.e 2,6,5,2 for 2^6 * 5^2  = 1600
   fn factor(&self)-> Vec<Self> where Self: Sized;
   /// Binary gcd, Computes the greatest common divisor of two numbers
   fn euclid_gcd(&self, other: &Self) -> Self;
   /// Counts the number of coprimes from 0 to self
   fn euler_totient(&self) -> Self;
     /// Computes self^p mod n
   fn mod_pow(&self, pow: &Self, n: &Self) -> Self;
   /// Computes the remainder of self*other mod n
   fn mul_mod(&self, other: &Self, n: &Self) -> Self;
    /// Identical to mul_modexcept in the case of u128 (not currently implemented) which uses a more optimized squaring
   fn quadratic_residue(&self, n: &Self) -> Self;
   /// Determines of a number is k-free, i.e square-free, cube-free etc. . .
   fn k_free(&self, k: &Self) -> bool;
   //fn partition(&self) -> u128;
   /// Returns the integer radical of self
   fn radical(&self)->Self;
   
   /// Smoothness bound of number
   fn smooth(&self) -> Self;
   /// Checks if the smoothness bound is at least b
   fn is_smooth(&self, b: &Self) -> bool;
   
   
   
   /// Legendre symbol of a,p. Assumes that p is an odd prime
   fn legendre(&self, p: &Self) -> i8 ;   
   /// Legendre symbol of a,p. Verifies that p is an odd prime, returns None if not. 
   fn checked_legendre(&self, p: &Self) -> Option<i8> ;
   
   
   /// Jacobi symbol of self,p. Assumes that p is an odd prime
   fn jacobi(&self, p: &Self) -> i8 ;  
   /// Jacobi symbol of self,p. Verifies that p is an odd prime, returns None if not. 
   fn checked_jacobi(&self, p: &Self) -> Option<i8> ;
   /*
   
   /// Pi approximation
   fn pi(&self) -> f64;
   
   let x = self.ln();
    //(self/x )*(1.0 + 1.0/x  + 2.0/(x.ln()*x.ln()))
       
   
   
   /// Louiville function
   fn louiville
   /// Mobius function
   fn mobius
   /// Mangoldt function 
   fn mangoldt
   /// Dedekind-psi
   fn dedekind_psi
   /// Jordan totient
   fn jordan_totient
   
   
   
    /// Integer partition
   fn partition()
   */
}

/*
NTT research 

  https://stackoverflow.com/questions/18465326/fast-bignum-square-computation
  
  https://stackoverflow.com/questions/18577076/modular-arithmetics-and-ntt-finite-field-dft-optimizations
  https://stackoverflow.com/questions/10260927/translation-from-complex-fft-to-finite-field-fft/18547575#18547575
  
  https://dl.acm.org/doi/10.1145/1185448.1185568

*/

