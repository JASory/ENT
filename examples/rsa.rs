use number_theory::Mpz;
use number_theory::NumberTheory;

/*

 Naive implementation of RSA cryptosystem in k-bits, works for k greater than 256 and

 This is not a cryptographic library and this implementation is not secure.

 Cryptographic libraries require special consideration for each operations and protections against variable manipulation

*/

// Generate safe prime
fn gen_safe(k: u32) -> Mpz {
    // Generate safe primes to minimize the factors of p-1 using standard primes are much faster

    loop {
        let p = Mpz::prime_gen(k);
        if p.congruence_u64(3, 2) {
            match p.is_sophie() {
                Some(x) => return x,
                None => (),
            }
        }
    }
}

fn generate_keys(k: u32) -> (Mpz, Mpz, Mpz) {
    // private, public, exponent
    let mut lo = k >> 1;
    let hi = lo;
    if k & 1 == 0 {
        // Ensures delta of atleast 2^(k-1) to prevent Fermat and Euler factorization
        lo -= 1;
    }
    let p = gen_safe(hi);
    //  println!("p = {}", p.to_string());
    let q = gen_safe(lo);
    //println!("q ={}", q.to_string());
    let phi = p
        .ref_subtraction(&Mpz::one())
        .ref_product(&q.ref_subtraction(&Mpz::one()));
    let n = p.ref_product(&q);
    assert_eq!(n.bit_length(), k as u64); // checks that the key is exactly k-bits

    let e = Mpz::from_u64(u64::prime_gen(64)); // generates exponent between 2^63;2^64
    let mut inv = e.extended_gcd(&phi).1;
    if !inv.is_positive() {
        inv = phi.ref_addition(&inv);
    }

    return (inv, n, e);
}

fn main() {
    const STRENGTH: u32 = 1024;

    let secret = Mpz::from_u64(0x29A822A1AA29229 << 5);
    let start = std::time::Instant::now();
    let (private, public, exponent) = generate_keys(STRENGTH);
    let stop = start.elapsed();
    let encrypted = secret.exp_residue(&exponent, &public);

    let decrypted = encrypted.exp_residue(&private, &public);
    println!("{:?} for a key of {}-bits \n ", stop, STRENGTH);
    println!(
        "{}",
        String::from_utf8(decrypted.to_u64().unwrap().to_be_bytes().to_vec()).unwrap()
    );
}
