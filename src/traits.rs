/*
  Traits

*/

pub trait NumberTheory{
   
   ///Returns optimally fast primality check for each interval

   fn is_prime(&self) -> bool;
   
   /// Factorizes into a vector of the form factor, power, factor, power . . . i.e 2,6,5,2 for 2^6 * 5^2  = 1600
   fn factor(&self)-> Vec<Self> where Self: Sized;
   /// Binary gcd, Computes the greatest common divisor of two numbers
   fn gcd(&self, other: &Self) -> Self;
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
   
   
}


