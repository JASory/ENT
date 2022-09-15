use number_theory::NumberTheory;

#[test]
fn product_residue() {
    assert_eq!(11u8.product_residue(&19, &0), 209);
    for i in 1..256u64 {
        assert_eq!(11u64.product_residue(&19, &i), 11 * 19 % i);
    }

    assert_eq!(17u16.product_residue(&19, &0), 323);
    for i in 1..256u16 {
        assert_eq!(17u16.product_residue(&19, &i), 17 * 19 % i);
    }

    assert_eq!(17u32.product_residue(&19, &0), 323);
    for i in 1..256u32 {
        assert_eq!(17u32.product_residue(&19, &i), 17 * 19 % i);
    }

    assert_eq!(17u64.product_residue(&19, &0), 323);
    for i in 1..256u64 {
        assert_eq!(17u64.product_residue(&19, &i), 17 * 19 % i);
    }
}

#[test]
fn exp_residue() {
    assert_eq!(0u64.checked_exp_residue(&0, &0).unwrap(), 1);
    assert_eq!(5u64.checked_exp_residue(&0, &0).unwrap(), 1);
    assert_eq!(5u64.checked_exp_residue(&1, &0).unwrap(), 5);
    assert_eq!(5u64.checked_exp_residue(&5, &0).unwrap(), 3125);
    assert_eq!(5u64.checked_exp_residue(&5, &5).unwrap(), 0);
    assert_eq!(5u64.checked_exp_residue(&5, &7).unwrap(), 3);
    assert_eq!((-5i64).checked_exp_residue(&5, &7).unwrap(), 4);
    assert_eq!((5i64).checked_exp_residue(&5, &7).unwrap(), 3);

    let add = (5i64).checked_exp_residue(&5, &7).unwrap();
    let inv = (-5i64).checked_exp_residue(&5, &7).unwrap();
    assert_eq!((add + inv) % 7, 0);

    let mul = (5i64).checked_exp_residue(&5, &7).unwrap();
    let inv = (5i64).checked_exp_residue(&-5, &7).unwrap();

    assert_eq!(mul * inv % 7, 1);
    // Unable to compute as inverse does not exist
    assert_eq!((5i64).checked_exp_residue(&-5, &35), None);
    // Unable to compute as the result exceeds the datatype
    assert_eq!((5i64).checked_exp_residue(&28, &0), None);
}
