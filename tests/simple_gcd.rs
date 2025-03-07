use number_theory::{Mpz,NumberTheory};

#[test]
fn gcds() {
    assert_eq!(35u8.gcd(63), 7);
    assert_eq!(35u16.gcd(63), 7);
    assert_eq!(35u32.gcd(63), 7);
    assert_eq!(35u64.gcd(63), 7);
    assert_eq!(35u128.gcd(63), 7);
    assert_eq!(
        Mpz::from(35u64)
            .gcd(Mpz::from(63u64)),
        Mpz::from(7u64)
    );

    let (gcd, bz1, bz2) = Mpz::from(36u64).extended_gcd(Mpz::from(81u64));
    assert_eq!(
        bz1.ref_product(&Mpz::from(36u64))
            .ref_addition(&bz2.ref_product(&Mpz::from(81u64))),
        gcd
    );
}
