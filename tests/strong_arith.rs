/*

  Randomly generated inputs to search for potential errors

  10-million integers in the interval -2^704-1 ; 2^704-1 are checked in the following operations takes approximately 5E+3 s

   Multiplication/Euclidean Division
   Logarithm
   Exponentiation/Nth-roots

*/
use number_theory::Mpz;
use number_theory::NumberTheory;

#[ignore]
#[test]
fn randarith() {
    const MAX_LOG_ERROR: f64 = 0.000001; // Maximum error allowed to be encurred by the product operation

    const ITERATIONS: usize = 10_000_000;

    const WORD_LENGTH: u64 = 11;

    let rand_gen = || -> Mpz {
        let mut z = Mpz::rand(WORD_LENGTH as usize); // Randomly generates integers in the range -2^704-1;2^704-1
        if u64::rng() % 2 == 0 {
            z.neg();
        }
        z
    };

    for _ in 0..ITERATIONS {
        let lhs = rand_gen(); // initializes first integer
        let rhs = rand_gen(); // initializes second integer
        let third = rand_gen(); // initializes third integer

        let lhsln = lhs.abs().ln(); // logarithms only apply to positive numbers
        let rhsln = rhs.abs().ln();
        let prodln_inf = lhsln + rhsln - MAX_LOG_ERROR; // Minimum bound allowed for the computed logarithm to take
        let prodln_sup = lhsln + rhsln + MAX_LOG_ERROR; // Maximum bound allowed for the computed logarithm to take

        assert_eq!(lhs, Mpz::from_string(&lhs.to_string()).unwrap()); // Checks that to_string and from_string are inverses
        assert_eq!(rhs, Mpz::from_string(&rhs.to_string()).unwrap());
        assert_eq!(third, Mpz::from_string(&third.to_string()).unwrap());

        let mut prod = lhs.ref_product(&rhs);
        let prodln = prod.abs().ln();

        assert!(prodln < prodln_sup && prodln > prodln_inf); // checks that the natural logarithm and the prod function behave roughly as expected

        prod.mut_addition(third);
        let (quo, rem) = prod.euclidean_div(&rhs);
        assert_eq!(quo.ref_product(&rhs).ref_addition(&rem), prod);

        assert_eq!(lhs.ref_product(&lhs).abs().sqrt().0, lhs.abs()); // checks the sqrt function, only applies to positive integers

        let randexp = (u64::rng() % 10) + 2;
        let estpowln = lhsln * (randexp as f64);
        let pow = lhs.pow(randexp);
        let powln = pow.abs().ln();

        let pow_inf = estpowln - MAX_LOG_ERROR * (randexp as f64);
        let pow_sup = estpowln + MAX_LOG_ERROR * (randexp as f64);

        assert!(powln < pow_sup && powln > pow_inf);

        assert_eq!(lhs.abs(), pow.abs().nth_root(&Mpz::from_u64(randexp)).0)
    }
}
