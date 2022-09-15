/*

  List of computationally intensive tests that demonstrate the claims made by the author.


  Full computation of each individual test may take hours so it is imperative that it is run as cargo test --release and during downtime.

  These tests are primarily to ensure confidence in program correctness, especially in light of the errors presented by other number-theory libraries in Rust.

  The assumptions that are made to claim correctness are listed next to the test. If the assumptions are wrong then the proof of correctness is also wrong

  Assumption abbreviations :  P = mathematically proven, V = Verified by at least two parties including at least partially by the author,
   UV = Computed by the author, but unverified by external parties,
   SP = statistically probable, applies to analyzing external software that corroborates well with proven bounds

*/

use number_theory::NumberTheory;

/*

   Description                      Assumptions                                                   Approximate Time
_______________________________________________________________________________________________________________
 Exhaustive proof of            Correctness of Erastothenes sieve in computing pi(n) (P)              5E+3 s
 deterministic primality  	Correctness of strong fermat test (P)
 for integers under 		Correctness of Feitsma-Galway pseudoprime table (V)
 2^64 + 2^44                    Correctness of J.A Sory's reduction and extension of table (UV)
                                Correctness of Kim Walisch's implementation of Gourdon's pi(n) (SP)
                                Correctness of Damgard et .al in pseudoprime density (P)


*/

/*
   Description                      Assumptions                                                   Approximate Time
___________________________________________________________________________________________________________________________

Disproof of Polya's                 Polya's conjecture is false (P)                              Factorization is currently too slow
conjecture and bound                Louiville function implementation is consistent
checking of  values of              regardless of the size of input (P)
Louiville function                 Factorization algorithm is correct (P)

*/

/*
fn liouville(){

let lvsum = let lsum = |x : u32| {
  let mut sum = -1i8;
   for j in 2..x+1{
    sum+=j.louiville();
   }
   return sum
   };

}
*/
/*
  Description                      Assumptions                                                   Approximate Time
___________________________________________________________________________________________________________________________

Proof of GCD correctness           GCD implementation is consistent


*/
#[ignore]
#[test]
fn strong_gcd() {
    let _euler_totient = 0u64;
    /*
    for i in 0..u32::MAX {
        if i.euclid_gcd(&u32::MAX) == 1 {
            euler_totient += 1;
        }
    }
    assert_eq!(2147483648, euler_totient);
    */
    for _ in 0..1_000_000_000 {
        // One billion iterations

        let au64 = u64::rng();
        let bu64 = u64::rng();

        let (gcd, a_inv, b_inv) = au64.extended_gcd(&bu64);

        if gcd == au64 || gcd == bu64 {
            // N extended gcd check
            assert_eq!((a_inv, b_inv), (0, 1))
        } else {
            assert_eq!(au64.product_residue(&a_inv, &bu64), gcd);
            assert_eq!(bu64.product_residue(&b_inv, &au64), gcd);
        }

        // Z extended gcd check

        let ai64 = i64::rng();
        let bi64 = i64::rng();

        let (gcd, x, y) = ai64.extended_gcd(&bi64);

        assert_eq!(ai64 * x + bi64 * y, gcd)
    }
}
/*

 Description                      Assumptions                                                   Approximate Time
___________________________________________________________________________________________________________________________

Proof of Extended Euclidean      gcd(a,b) = a*x + b*y (P)
correctness                      EEA implementation is consistent



 Description                      Assumptions                                                   Approximate Time
___________________________________________________________________________________________________________________________
*/

/*
 Description                      Assumptions                                                   Approximate Time
___________________________________________________________________________________________________________________________

*/
#[ignore]
#[test]
fn factor_test() {
    // extremely slow, factorization must be improved
    let factorprod = |x: Vec<u64>| -> u64 {
        let mut prod = 1;

        for i in 0..x.len() / 2 {
            prod *= x[i * 2].pow(x[i * 2 + 1] as u32)
        }
        return prod;
    };

    for i in 1u64..u64::MAX >> 44 {
        assert_eq!(i, factorprod(i.factor()))
    }

    for _i in 0..1000 {
        // 60s critical to improve
        let rand = u64::rng();
        assert_eq!(rand, factorprod(rand.factor()))
    }
}
