# ENT

Elementary Number Theory for Integers in Rust

Algebraic definitions of primality and factorization are used, permitting checks like -127.is_prime() to return true and unique factorizations to be considered unsigned.Published as number-theory on crates.io

Currently implements these functions

- Primality
- Factorization
- GCD, Extended GCD
- Carmichael,Euler & Jordan totients 
- Dedekind psi
- Liouville, and Mobius function
- Prime-counting function/nth-prime, and prime lists
- Integer sqrt/nth root
- Integer radical
- K-free
- Quadratic and Exponential residues
- Legendre,Jacobi and Kronecker symbols
- Smoothness checks

 Additionally this library has an implementation of the previous NT functions for arbitrary-precision integers, plus some elementary arithmetic.
 Multiplication utilizes Karatsuba algorithm, otherwise all other arithmetic can be assumed to be naive. 
 
 - Addition/subtraction
 - Multiplication 
 - Euclidean Division 
 - Conversion to and from radix-10 string
 - Successor function (+1)
 - SIRP-factorials {generalization of factorials}
 - Conditional Interval Product (CIP factorial)
 - Sqrt/nth root
 - Exponentiation
 - Logarithms
 - Probable pseudoprime construction

Usage is fairly simple
 ```rust
 // include NT functions
 use number_theory::NumberTheory;
 // include arbitrary-precision arithmetic
 use number_theory::Mpz;
   // Sign, generally unnecessary for ENT
 //use number_theory:Sign; 
 let mersenne = Mpz::from_string("-127").unwrap(); 
 assert_eq!(mersenne.is_prime(), true);
  // Or for a more complex example
  
  // check if x mod 1 === 0, trivial closure
  let modulo = |x: &u64| {if x%1==0{return true} return false};
  
   //Generate  872 factorial, using the above trivial function
   // this can be just as easily reproduced as Mpz::sirp(1,872,1,0);
  let mut factorial = Mpz::cip(1, 872,modulo);
  
  // Successor function, increment by one
  factorial.successor();
  
  // 872! + 1 is in-fact a factorial prime
  assert_eq!(factorial.is_prime(),true)
 ```
