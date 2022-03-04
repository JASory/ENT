# ENT

Elementary Number Theory for Integers in Rust

Currently the fastest library for factorization and primality checking in the interval 0;2^64 that is available in Crates.io and possibly in the entire Rust-lang ecosystem. Algebraic definitions of primality and factorization are used, permitting checks like -127.is_prime() to return true and unique factorizations to be considered unsigned.



Currently implements these functions

- Primality
- Factorization
- GCD
- Euler totient
- Integer radical (not sqrt)
- K-free
- Modular exponentiation and quadratic residues
- Legendre symbol
- Jacobi symbol
- Smoothness

 Additionally this library has an implementation of the previous NT functions for arbitrary-precision integers, plus some elementary arithmetic. Multiplication utilizes Karatsuba algorithm, otherwise all other arithmetic can be assumed to be naive. 
 
 - Addition/subtraction
 - Multiplication
 - Euclidean Division
 - Extended Euclidean (EEA)
 - Conversion to and from radix-10 string
 - Successor function (+1)
 - SIRP-factorials
 - Conditional Interval Product (CIP factorial)
 - Sqrt/nth root
 - Exponentiation
 - Logarithms

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
  let mut factorial = Mpz::cip(1, 872,modulo);
  
  // Successor function, increment by one
  factorial.successor();
  
  // 872! + 1 is in-fact a factorial prime
  assert_eq!(factorial.is_prime(),true)
 ```

[![GitHub stats](https://github-readme-stats.vercel.app/api?username=JASory)](https://github.com/JASory/github-readme-stats)
