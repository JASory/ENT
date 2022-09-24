/*!
A  library for elementary number theory and arbitrary-precision arithmetic.Intended to be simple,  versatile, performant and above all else correct.

 Unlike most other number-theoretic libraries number-theory is more generalized
operating over quotient rings (also called residue-classes in the methods, i.e exp_residue for the residue class of pow(x,y) computed in Z/nZ), rather than using modular arithmetic definitions, this permits more general results which in some cases is different than the typical modular arithmetic definition. We provide most of them below

* Z/0Z = Z , consequently arithmetic performed over it is equivalent to arithmetic in Z, hence 4.quadratic_residue(0) == 8
* Z/-nZ = Z/nZ  and therefore  4.quadratic_residue(-3) = 4.quadratic_residue(3).
* Negative numbers are evaluated as additive inverses in the ring Z/nZ for simplicity (except in the case of Z/0Z), so (-4).exp_residue(&3,5) == 1 not -4 although they are equivalent

Features

* 14 distinct non-trivial number-theorectic functions implemented by all builtin integer types and arbitrary precision integer
* Z/nZ Ring arithmetic, frequently using optimized algorithms.
* Extremely fast and rigourously proven primality checking for integers under 2^64+2^45.

Considerations for usage

 * Smallest datatype is usually most efficient, both due to the processor and the fact that many of the functions are better optimized for small values.
 * Functions that have a "checked" counterpart perform zero correctness checks, **do not panic**, and rely on the user guaranteeing that the input will be valid.
 This is for efficiency purposes, if you cannot guarantee that an input will meet the guidelines then use the "checked" variant. The inputs that result in incorrect  outputs are given under "Failure" for each applicable function. This conditions directly correspond to the conditions that result in None output in the checked variant. 
 * Number-theory is not a cryptography library and is not a substitute for rigourous securely implemented algorithms that consider adversarial attacks, branch-prediction etc. However, number-theory may exceed other "cryptographic" libraries in speed and correctness, possibly even under adversarial conditions.

Why you might not want to use this

 * API will break without warning. Number-theory is not even close to a stable release it is simply published to make it available given the fact that it surpasses many similar libraries. (Indeed number-theory was originally a subfile to a still private computer algebra project, and consequently it's functionality strongly reflects this).
 * Functions may be silently broken. Outside of is_prime no functions have a rigourous proof of correctness, (although few other libraries do either, indeed "silently broken" is the norm in mathematical software developed by non-specialists).
 * It doesn't use dependencies or follow any external api system (i.e num-traits). While deliberate to maintain simplicity and greater control of the software on the developer side, users may find it less versatile.
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

pub use crate::arithmetic::mpz::Mpz;
pub use crate::arithmetic::sign::Sign;
pub use crate::traits::NumberTheory;
