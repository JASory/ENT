 use number_theory::Mpz;
 use number_theory::NumberTheory;
 use number_theory::Sign; 
   
   //#[ignore]
   #[test]
    fn primality(){
    
    
    assert_eq!(0u8.is_prime(),false);
    assert_eq!(2u8.is_prime(),true);
    
     assert_eq!(0i8.is_prime(),false);
    assert_eq!(2i8.is_prime(),true);
    

  let mut mersenne = Mpz::from_u64(1).shl(82589933);  //1290467  44
     mersenne.mut_subtraction(Mpz::one());
     
     assert_eq!(mersenne.is_prime(),true)
    
    }
    
    //#[ignore]
    #[test]
    fn euclidean_div(){
    assert_eq!(17i8.euclidean_div(&5i8).0 * 5i8 + (17i8.euclidean_div(&5i8).1), 17i8);
    assert_eq!(17i8.euclidean_div(&-5i8).0 * -5i8 + (17i8.euclidean_div(&-5i8).1), 17i8);
    assert_eq!(-17i8.euclidean_div(&-5i8).0 * -5i8 + (-17i8.euclidean_div(&-5i8).1), -17i8);
    assert_eq!(-17i8.euclidean_div(&5i8).0 * 5i8 + (-17i8.euclidean_div(&5i8).1), -17i8);
    
    assert_eq!(17i16.euclidean_div(&5i16).0 * 5i16 + (17i16.euclidean_div(&5i16).1), 17i16);
    assert_eq!(17i16.euclidean_div(&-5i16).0 * -5i16 + (17i16.euclidean_div(&-5i16).1), 17i16);
    assert_eq!(-17i16.euclidean_div(&-5i16).0 * -5i16 + (-17i16.euclidean_div(&-5i16).1), -17i16);
    assert_eq!(-17i16.euclidean_div(&5i16).0 * 5i16 + (-17i16.euclidean_div(&5i16).1), -17i16);
    
    assert_eq!(17i32.euclidean_div(&5i32).0 * 5i32 + (17i32.euclidean_div(&5i32).1), 17i32);
    assert_eq!(17i32.euclidean_div(&-5i32).0 * -5i32 + (17i32.euclidean_div(&-5i32).1), 17i32);
    assert_eq!(-17i32.euclidean_div(&-5i32).0 * -5i32 + (-17i32.euclidean_div(&-5i32).1), -17i32);
    assert_eq!(-17i32.euclidean_div(&5i32).0 * 5i32 + (-17i32.euclidean_div(&5i32).1), -17i32);
    
    assert_eq!(17i64.euclidean_div(&5i64).0 * 5i64 + (17i64.euclidean_div(&5i64).1), 17i64);
    assert_eq!(17i64.euclidean_div(&-5i64).0 * -5i64 + (17i64.euclidean_div(&-5i64).1), 17i64);
    assert_eq!(-17i64.euclidean_div(&-5i64).0 * -5i64 + (-17i64.euclidean_div(&-5i64).1), -17i64);
    assert_eq!(-17i64.euclidean_div(&5i64).0 * 5i64 + (-17i64.euclidean_div(&5i64).1), -17i64);
    
    assert_eq!(17i128.euclidean_div(&5i128).0 * 5i128 + (17i128.euclidean_div(&5i128).1), 17i128);
    assert_eq!(17i128.euclidean_div(&-5i128).0 * -5i128 + (17i128.euclidean_div(&-5i128).1), 17i128);
    assert_eq!(-17i128.euclidean_div(&-5i128).0 * -5i128 + (-17i128.euclidean_div(&-5i128).1), -17i128);
    assert_eq!(-17i128.euclidean_div(&5i128).0 * 5i128 + (-17i128.euclidean_div(&5i128).1), -17i128);
    
    
    assert_eq!(17i16.euclidean_div(&5i16).0 * 5i16 + (17i16.euclidean_div(&5i16).1), 17i16);
    assert_eq!(17i16.euclidean_div(&-5i16).0 * -5i16 + (17i16.euclidean_div(&-5i16).1), 17i16);
    assert_eq!(-17i16.euclidean_div(&-5i16).0 * -5i16 + (-17i16.euclidean_div(&-5i16).1), -17i16);
    assert_eq!(-17i16.euclidean_div(&5i16).0 * 5i16 + (-17i16.euclidean_div(&5i16).1), -17i16);
    
    let (quo, rem) = Mpz::from_i64(17).euclidean_div(&Mpz::from_i64(5));
    assert_eq!( (quo.ref_product(&Mpz::from_i64(5))).ref_addition(&rem), Mpz::from_i64(17));
    
    let (quo, rem) = Mpz::from_i64(-17).euclidean_div(&Mpz::from_i64(5));
    assert_eq!( (quo.ref_product(&Mpz::from_i64(5))).ref_addition(&rem), Mpz::from_i64(-17));
    
    
     let (quo, rem) = Mpz::from_i64(-17).euclidean_div(&Mpz::from_i64(-5));
    assert_eq!( (quo.ref_product(&Mpz::from_i64(-5))).ref_addition(&rem), Mpz::from_i64(-17));
    
     let (quo, rem) = Mpz::from_i64(17).euclidean_div(&Mpz::from_i64(-5));
    assert_eq!( (quo.ref_product(&Mpz::from_i64(-5))).ref_addition(&rem), Mpz::from_i64(17));
    // i128                                              
    let (quo, rem) = Mpz::from_i128(17).euclidean_div(&Mpz::from_i128(5));
    assert_eq!( (quo.ref_product(&Mpz::from_i128(5))).ref_addition(&rem), Mpz::from_i128(17));
    
    let (quo, rem) = Mpz::from_i128(-170000000000000000000).euclidean_div(&Mpz::from_i128(50000000));
    assert_eq!( (quo.ref_product(&Mpz::from_i128(50000000))).ref_addition(&rem), Mpz::from_i128(-170000000000000000000));
    
    let (quo,rem) = Mpz::from_i128(-5i128).euclidean_div(&Mpz::from_i128(1i128));
    assert_eq!(quo.to_string(),"-5");
    assert_eq!(rem.to_string(),"0");
    
    let (quo,rem) = Mpz::from_string("-2500000000000000000000000000000000000000000000000000000000000000000000000000000000").unwrap().euclidean_div(&Mpz::from_string("50000000000000000000000000000000000000000").unwrap());
    assert_eq!(rem.to_string(),"0");
    
    let (quo,rem) = Mpz::from_string("50000000000000000000000000000000000000000").unwrap().euclidean_div(&Mpz::from_string("-2500000000000000000000000000000000000000000000000000000000000000000000000000000000").unwrap());
    println!("{:?}",quo);
    assert_eq!(quo,Mpz::unchecked_new(Sign::Positive,vec![0]));
    }
    
   
   
