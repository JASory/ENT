# ENT

Elementary Number Theory for Integers in Rust

Currently the fastest library for factorization and primality checking in the interval 0;2^64 that is available in Crates.io and possibly in the entire Rust-lang ecosystem. Algebraic definitions of primality and factorization are used, permitting checks like -127.is_prime() to return true and unique factorizations to be considered unsigned.



Currently implements these functions

- Primality
- Factorization
- Euler totient
- Integer radical (not sqrt)
- K-free
- Modular exponentiation and quadratic residues
- Legendre symbol
- Jacobi symbol

 Additionally this library has an implementation of the previous NT functions for arbitrary-precision integers, plus some elementary arithmetic. Multiplication utilizes Karatsuba algorithm, otherwise all other arithmetic can be assumed to be naive. 
 
 - Addition/subtraction
 - Multiplication
 - Euclidean Division
 - Conversion to and from radix-10 string
 - Successor function (+1)
 - SIRP-factorials
 - Sqrt/nth root
 - Exponentiation
 - Logarithms

Usage is fairly simple
 ```rust
 // include NT functions
 use number_theory::traits::NumberTheory;  
 // include arbitrary-precision arithmetic
 use number_theory::arithmetic::mpz::Mpz;
   // Sign, generally unnecessary for ENT
 //use number_theory::arithmetic::sign::Sign; 
 // unsigned from string, defaults to Sign::Positive
 let mersenne = Mpz::from_string("-127"); 
 assert_eq!(mersenne.is_prime(), true);
 ```
