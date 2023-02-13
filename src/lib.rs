/*!
A  library for elementary number theory and arbitrary-precision integer arithmetic.Intended to be simple,  versatile, performant and above all else
 correct.

 Unlike most other number-theoretic libraries number-theory is more generalized operating over quotient rings (also called residue-classes
 in the methods, i.e exp_residue for the residue class of pow(x,y) computed in Z/nZ), rather than using modular arithmetic definitions, this
 permits more general results which in some cases is different than the typical modular arithmetic definition. We provide most of them below

* Z/0Z = Z , consequently arithmetic performed over it is equivalent to arithmetic in Z, hence 4.quadratic_residue(0) == 16
* Z/-nZ = Z/nZ  and therefore  4.quadratic_residue(-3) = 4.quadratic_residue(3).
* Negative numbers are evaluated as additive inverses in the ring Z/nZ for simplicity (except in the case of Z/0Z), so (-4).exp_residue(&3,5) == 1
 not -4 although they are equivalent
* Zero is considered to have infinite factors, 1 is considered to have zero prime factors. These conditions result in the need for checked variants
of many functions which are more useful in handling the 0 and 1 cases even if the unchecked variants are strictly correct. 

Features

* 14 distinct non-trivial number-theorectic functions implemented by all builtin integer types and arbitrary precision integer (some of the
  provided functions are indeed relatively trivial or strongly derived from others and are so not counted here)
* Z/nZ Ring arithmetic, frequently using optimized algorithms.
* Extremely fast and rigourously proven primality checking for integers under 2^64+2^49.

Considerations for usage

 * Smallest datatype is usually most efficient, both due to the processor and the fact that many of the functions are better optimized for small
  values.
 * Functions that have a "checked" counterpart perform zero correctness checks, **do not panic**, and rely on the user guaranteeing that the input
  will be valid.This is for efficiency purposes, if you cannot guarantee that an input will meet the guidelines then use the "checked" variant. 
  The inputs that result in incorrect  outputs are given under "Failure" for each applicable function. This conditions directly correspond to the
  conditions that result in  output in the checked variant. 
 * As Number-theory functions are often more generalized than most other libraries this means than some common functionality may be hidden in
  related  functions. Examples include multiplicative inverses being provided by x.exp_residue(-1,n) or more directly by the finite ring variant
  of euclid_gcd, and detecting prime powers using the mangoldt function or max_root
  
 * Number-theory is not a cryptography library and is not a substitute for rigourous securely implemented algorithms that consider adversarial
  attacks, branch-prediction etc. However, number-theory may exceed other "cryptographic" libraries in speed and correctness, possibly even under
  adversarial conditions. Additionally it provides functions that may be useful in testing cryptography libraries like 
  [psp](struct.Mpz.html#method.psp) which provides a composite with optimal probability of passing a Miller-Rabin test. 
  
Rules of Checked Functions and NTResult

 * Checked functions exist to provide a function that returns a correct evaluation for any input at the cost of greater overhead 
 * Checked Functions always return an NTResult enum
 * Checked and non-checked variants of functions exist if the function: 1. Cannot correctly evaluate all inputs 2. Has low complexity
 * "Non-checked" functions, without a checked counterpart, may return an NTResult if they: 1. Can error  2. Have high complexity. The majority of
   these functions rely on costly factorization and evaluate as infinite for 0. 

Why you might not want to use this

 * API will break without warning. Number-theory is not even close to a stable release it is simply published to make it available given the fact
  that it surpasses many similar libraries.
 
 * Functions may be silently broken. Outside of is_prime no functions have a rigourous proof of correctness, (although few other libraries do
  either, indeed "silently broken" is the norm in mathematical software developed by non-specialists).
 * It doesn't use dependencies or follow any external api system (i.e num-traits). While deliberate to maintain simplicity and greater control of
  the software on the developer side, users may find it less versatile.
 * Only works on x86-64 platforms

 Future Work
 * Subquadratic division
 * Gourdon's pi
 * Switch from N-1 primality proofing to Elliptic Curves
 * Implement GNFS
 * Prove correctness for all functions
 
 Only once the last is achieved will number-theory be stabilized. 



*/

/* lints to use
missing_docs, trivial_casts, trivial_numeric_casts,
unused_import_braces, unused_qualifications
*/
#![forbid(warnings,unused_parens, missing_docs,trivial_casts, trivial_numeric_casts,arithmetic_overflow)]

pub(crate) mod arithmetic;
mod data;
pub(crate) mod montgomery;
mod primitive;
mod sieve;
mod traits;
mod result;

pub use crate::arithmetic::mpz::Mpz;
pub use crate::arithmetic::sign::Sign;
pub use crate::traits::NumberTheory;
pub use crate::result::NTResult;
