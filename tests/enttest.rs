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
    
    
    
    //assert_eq!()
    assert_eq!(Mpz::u_new(vec![
14125993214864411435, 17627421051211808665, 8803062915072472141, 16592776965248599063,
2400842081300231610, 15499637576960551660, 13635884493322354475, 11162303294968070407, 
16638949889857351534, 4755447895431218293, 3008559606663544016, 6962752268970074012, 
15476370627519396605, 7818823728409314089, 7726208862163452462, 2793947945622957183, 
9898516143788448943, 2723958373070846137, 15945921893246334102, 9257474160255769327, 138699416178]).probable_prime(), false);

  let mut mersenne = Mpz::from_u64(1).shl(82589933);  //1290467  44
     mersenne.mut_subtraction(Mpz::one());
     
     assert_eq!(mersenne.probable_prime(),true)
    
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
    
   
   
