use number_theory::Mpz;
use number_theory::NumberTheory;
use number_theory::Sign;

#[test]
fn primality() {
    // proves that integers below 2^16 are correctly evaluated
    assert_eq!(0u8.is_prime(), false);
    assert_eq!(2u8.is_prime(), true);

    assert_eq!(0i8.is_prime(), false);
    assert_eq!(2i8.is_prime(), true);

    let mut mersenne = Mpz::from_u64(1).shl(82589933); //1290467  44
    mersenne.mut_subtraction(Mpz::one());

    assert_eq!(mersenne.is_prime(), true);

    // prime ranges, for checks beyond this value see det_primality in strong_prime which proves up to 2^64+2^42

    let mut count = 0u32;

    for i in 0..u8::MAX {
        if i.is_prime() {
            count += 1
        }
    }
    assert_eq!(count, 54u32);

    for i in 255..u16::MAX {
        if i.is_prime() {
            count += 1
        }
    }

    assert_eq!(count, 6542u32);

    for i in -65535i32..0i32 {
        // Demonstrate that negative primes work correctly as well
        if i.is_prime() {
            count -= 1;
        }
    }
    assert_eq!(count, 0);
    let sevenbit = i8::prime_gen(7) as u8;
    assert!((sevenbit < 127) && (sevenbit > 64));
    
    let thirtybit = u32::prime_gen(30);
    assert!((thirtybit < 1073741824) && (thirtybit > 536870912))
}

#[test]
fn euclidean_div() {
    assert_eq!(
        17i8.euclidean_div(&5i8).0 * 5i8 + (17i8.euclidean_div(&5i8).1),
        17i8
    );
    assert_eq!(
        17i8.euclidean_div(&-5i8).0 * -5i8 + (17i8.euclidean_div(&-5i8).1),
        17i8
    );
    assert_eq!(
        -17i8.euclidean_div(&-5i8).0 * -5i8 + (-17i8.euclidean_div(&-5i8).1),
        -17i8
    );
    assert_eq!(
        -17i8.euclidean_div(&5i8).0 * 5i8 + (-17i8.euclidean_div(&5i8).1),
        -17i8
    );

    assert_eq!(
        17i16.euclidean_div(&5i16).0 * 5i16 + (17i16.euclidean_div(&5i16).1),
        17i16
    );
    assert_eq!(
        17i16.euclidean_div(&-5i16).0 * -5i16 + (17i16.euclidean_div(&-5i16).1),
        17i16
    );
    assert_eq!(
        -17i16.euclidean_div(&-5i16).0 * -5i16 + (-17i16.euclidean_div(&-5i16).1),
        -17i16
    );
    assert_eq!(
        -17i16.euclidean_div(&5i16).0 * 5i16 + (-17i16.euclidean_div(&5i16).1),
        -17i16
    );

    assert_eq!(
        17i32.euclidean_div(&5i32).0 * 5i32 + (17i32.euclidean_div(&5i32).1),
        17i32
    );
    assert_eq!(
        17i32.euclidean_div(&-5i32).0 * -5i32 + (17i32.euclidean_div(&-5i32).1),
        17i32
    );
    assert_eq!(
        -17i32.euclidean_div(&-5i32).0 * -5i32 + (-17i32.euclidean_div(&-5i32).1),
        -17i32
    );
    assert_eq!(
        -17i32.euclidean_div(&5i32).0 * 5i32 + (-17i32.euclidean_div(&5i32).1),
        -17i32
    );

    assert_eq!(
        17i64.euclidean_div(&5i64).0 * 5i64 + (17i64.euclidean_div(&5i64).1),
        17i64
    );
    assert_eq!(
        17i64.euclidean_div(&-5i64).0 * -5i64 + (17i64.euclidean_div(&-5i64).1),
        17i64
    );
    assert_eq!(
        -17i64.euclidean_div(&-5i64).0 * -5i64 + (-17i64.euclidean_div(&-5i64).1),
        -17i64
    );
    assert_eq!(
        -17i64.euclidean_div(&5i64).0 * 5i64 + (-17i64.euclidean_div(&5i64).1),
        -17i64
    );

    assert_eq!(
        17i128.euclidean_div(&5i128).0 * 5i128 + (17i128.euclidean_div(&5i128).1),
        17i128
    );
    assert_eq!(
        17i128.euclidean_div(&-5i128).0 * -5i128 + (17i128.euclidean_div(&-5i128).1),
        17i128
    );
    assert_eq!(
        -17i128.euclidean_div(&-5i128).0 * -5i128 + (-17i128.euclidean_div(&-5i128).1),
        -17i128
    );
    assert_eq!(
        -17i128.euclidean_div(&5i128).0 * 5i128 + (-17i128.euclidean_div(&5i128).1),
        -17i128
    );

    assert_eq!(
        17i16.euclidean_div(&5i16).0 * 5i16 + (17i16.euclidean_div(&5i16).1),
        17i16
    );
    assert_eq!(
        17i16.euclidean_div(&-5i16).0 * -5i16 + (17i16.euclidean_div(&-5i16).1),
        17i16
    );
    assert_eq!(
        -17i16.euclidean_div(&-5i16).0 * -5i16 + (-17i16.euclidean_div(&-5i16).1),
        -17i16
    );
    assert_eq!(
        -17i16.euclidean_div(&5i16).0 * 5i16 + (-17i16.euclidean_div(&5i16).1),
        -17i16
    );

    let (quo, rem) = Mpz::from_i64(17).euclidean_div(&Mpz::from_i64(5));
    assert_eq!(
        (quo.ref_product(&Mpz::from_i64(5))).ref_addition(&rem),
        Mpz::from_i64(17)
    );

    let (quo, rem) = Mpz::from_i64(-17).euclidean_div(&Mpz::from_i64(5));
    assert_eq!(
        (quo.ref_product(&Mpz::from_i64(5))).ref_addition(&rem),
        Mpz::from_i64(-17)
    );

    let (quo, rem) = Mpz::from_i64(-17).euclidean_div(&Mpz::from_i64(-5));
    assert_eq!(
        (quo.ref_product(&Mpz::from_i64(-5))).ref_addition(&rem),
        Mpz::from_i64(-17)
    );

    let (quo, rem) = Mpz::from_i64(17).euclidean_div(&Mpz::from_i64(-5));
    assert_eq!(
        (quo.ref_product(&Mpz::from_i64(-5))).ref_addition(&rem),
        Mpz::from_i64(17)
    );
    // i128
    let (quo, rem) = Mpz::from_i128(17).euclidean_div(&Mpz::from_i128(5));
    assert_eq!(
        (quo.ref_product(&Mpz::from_i128(5))).ref_addition(&rem),
        Mpz::from_i128(17)
    );

    let (quo, rem) =
        Mpz::from_i128(-170000000000000000000).euclidean_div(&Mpz::from_i128(50000000));
    assert_eq!(
        (quo.ref_product(&Mpz::from_i128(50000000))).ref_addition(&rem),
        Mpz::from_i128(-170000000000000000000)
    );

    let (quo, rem) = Mpz::from_i128(-5i128).euclidean_div(&Mpz::from_i128(1i128));
    assert_eq!(quo.to_string(), "-5");
    assert_eq!(rem.to_string(), "0");

    let (_quo, rem) = Mpz::from_string(
        "-2500000000000000000000000000000000000000000000000000000000000000000000000000000000",
    )
    .unwrap()
    .euclidean_div(&Mpz::from_string("50000000000000000000000000000000000000000").unwrap());
    assert_eq!(rem.to_string(), "0");

    let (quo,_rem) = Mpz::from_string("50000000000000000000000000000000000000000").unwrap().euclidean_div(&Mpz::from_string("-2500000000000000000000000000000000000000000000000000000000000000000000000000000000").unwrap());
    println!("{:?}", quo);
    assert_eq!(quo, Mpz::unchecked_new(Sign::Positive, vec![0]));
}

#[rustfmt::skip]
#[test]
fn factor() {
    assert_eq!(0xFFFFFFF600000019u64.factor(), vec![4294967291, 2]);
    assert_eq!(
        u64::MAX.factor(),
        vec![3, 1, 5, 1, 17, 1, 257, 1, 641, 1, 65537, 1, 6700417, 1]
    );
    assert_eq!(
        u128::MAX.factor(),
        vec![3,1,5,1,17,1,257,1,641,1,65537,1,274177,1,6700417,1,67280421310721,1]
    )
}

#[test]
fn k_free() {
    assert_eq!(108u8.k_free(&4), true);
    assert_eq!(9375u16.k_free(&5), false);

    let factorial = Mpz::sirp(1, 1000, 1, 0);

    for i in 2..15 {
        // Checks the power of the first 15 prime factors
        assert_eq!(factorial.k_free(&Mpz::from_u64(i)), false)
    }
}
