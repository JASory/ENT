use number_theory::Mpz;
use number_theory::NumberTheory;

#[test]
fn gcds() {
    assert_eq!(35u8.euclid_gcd(&63), 7);
    assert_eq!(35u16.euclid_gcd(&63), 7);
    assert_eq!(35u32.euclid_gcd(&63), 7);
    assert_eq!(35u64.euclid_gcd(&63), 7);
    assert_eq!(35u128.euclid_gcd(&63), 7);
    assert_eq!(
        Mpz::from_u64(35)
            .euclid_gcd(&Mpz::from_u64(63))
            .to_u64()
            .unwrap(),
        7
    );

    let (gcd, bz1, bz2) = Mpz::from_u64(36).extended_gcd(&Mpz::from_u64(81));
    assert_eq!(
        bz1.ref_product(&Mpz::from_u64(36))
            .ref_addition(&bz2.ref_product(&Mpz::from_u64(81))),
        gcd
    );
}
