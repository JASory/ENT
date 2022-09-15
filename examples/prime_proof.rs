use number_theory::NumberTheory;

/*

  Demonstrating how to verify primes from the prime_proof function

  The prime_proof function returns a claim, a witness and the factorization of N-1 (also known as a Pratt certificate)

  Verification requires that prime_proof's claim (the first boolean value) equals the verifier result

*/

// Input : an integer N and a Pratt Primality Certificate

fn verify(n: u64, p: Vec<u64>) -> bool {
    let witness = p[0];

    if witness.exp_residue(&(n - 1), &n) != 1 {
        // check a^n-1 mod n
        return false;
    }

    for i in 1..p.len() {
        let test = witness.exp_residue(&((n - 1) / p[i]), &n); // check a^(n-1)/p mod n if it equals 1 then the witness verifies that it is composite
        if test == 1 {
            return false;
        }
    }
    return true;
}

fn main() {
    let p = u64::rng();
    let (claim, certificate) = p.prime_proof();

    println!("\n \"{} is prime\" is a {} statement \n", p, claim);

    println!(
        "We have  a witness {} and the factors of {}-1 are as follows {:?} \n  ",
        certificate[0],
        p,
        &certificate[1..]
    );

    println!("By checking the cerificate we can see that the statement \"{} is prime \" is a {} statement  and that the certificate is {} \n ", 
  p, verify(p,certificate.clone()), verify(p,certificate)==claim)
}
