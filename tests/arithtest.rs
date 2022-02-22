use number_theory::Mpz;
use number_theory::Sign; 
 //#[ignore]
 #[test]
 
 fn arith(){
 
 assert_eq!(Mpz::zero().ref_addition(&Mpz::u_from_string("222222222222222222222").unwrap()).to_string(), "222222222222222222222");
 assert_eq!(Mpz::one().ref_addition(&Mpz::u_from_string("222222222222222222222").unwrap()).to_string(), "222222222222222222223");
 
 assert_eq!(Mpz::u_from_string("222222222222222222222").unwrap().ref_addition(&Mpz::u_from_string("222222222222222222223").unwrap()).to_string(), "444444444444444444445");
 assert_eq!(Mpz::u_from_string("222222222222222222222").unwrap().ref_addition(&
 Mpz::u_from_string("222222222222222222222222222222222222222222").unwrap()).to_string(), "222222222222222222222444444444444444444444");
 
 assert_eq!(Mpz::u_from_string("222222222222222222222222222222222222222222").unwrap().ref_addition(&
 Mpz::u_from_string("222222222222222222222").unwrap()).to_string(), "222222222222222222222444444444444444444444");
 
 assert_eq!(Mpz::from_string("-222222222222222222222222222222222222222222").unwrap().ref_addition(&
 Mpz::u_from_string("222222222222222222222").unwrap()).to_string(), "-222222222222222222222000000000000000000000");
 
 assert_eq!(Mpz::from_string("-222222222222222222222222222222222222222222").unwrap().ref_addition(&
 Mpz::from_string("-222222222222222222222").unwrap()).to_string(), "-222222222222222222222444444444444444444444");
 
 assert_eq!(Mpz::from_string("222222222222222222222222222222222222222222").unwrap().ref_addition(&
 Mpz::from_string("-222222222222222222222").unwrap()).to_string(), "222222222222222222222000000000000000000000");
 
  // Product                                                                                                   
 assert_eq!(Mpz::zero().ref_product(&Mpz::u_new(vec![1,2])),Mpz::zero());
 assert_eq!(Mpz::one().ref_product(&Mpz::u_new(vec![1,2])),Mpz::u_new(vec![1,2]));
 
 assert_eq!(Mpz::one().ref_product(&Mpz::u_new(vec![1,2])),Mpz::u_new(vec![1,2,0]));
  
 assert_eq!(Mpz::u_from_string("222222222222222222222").unwrap().ref_product(&Mpz::u_from_string("222222222222222222222").unwrap()).to_string(),
 "49382716049382716049283950617283950617284" ) ;
 
 assert_eq!(Mpz::from_string("-222222222222222222222").unwrap().ref_product(&Mpz::from_string("-222222222222222222222").unwrap()).to_string(),
 "49382716049382716049283950617283950617284") ;
 
 assert_eq!(Mpz::from_string("222222222222222222222").unwrap().ref_product(&Mpz::from_string("-222222222222222222222").unwrap()).to_string(),
 "-49382716049382716049283950617283950617284") ;
 
 assert_eq!(Mpz::from_string("-222222222222222222222").unwrap().ref_product(&Mpz::from_string("222222222222222222222").unwrap()).to_string(),
 "-49382716049382716049283950617283950617284") ;
 // exponentiation
 assert_eq!(Mpz::from_i64(33).pow(42).to_string(), "5992188611668647482797702362673386216879593110712943882639341889");
 assert_eq!(Mpz::from_i64(-33).pow(42).to_string(), "5992188611668647482797702362673386216879593110712943882639341889");
 assert_eq!(Mpz::from_i64(-33).pow(43).to_string(), "-197742224185065366932324177968221745157026572653527148127098282337");
 
 
 
 }
 
 // #[ignore]
 #[test]
 
 fn math(){
 
 assert_eq!(Mpz::from_u64(33).pow(43).sqrt().to_string(),"444682160857690992367099165368543"); // sqrt
 assert_eq!(Mpz::from_u64(33).pow(43).ln(),150.34982514305864);// natural logarithm
 assert_eq!(Mpz::from_u64(33).pow(43).log(2.718281828459045235),150.34982514305864);// natural logarithm
 
 //assert_eq!(Mpz::from_u64(2).shl(65536).iter_log(2.0),4);
 
 
 }
 //#[ignore]
 #[test] 
 
 fn bitwise(){
 // Shift-right ops
 assert_eq!(Mpz::from_u64(1).shl(0).to_u64().unwrap(),1);
 assert_eq!(Mpz::from_u64(1).shl(1).to_u64().unwrap(),2);
 assert_eq!(Mpz::from_u64(1).shl(65).to_u128().unwrap(),36893488147419103232);
 assert_eq!(Mpz::u_from_string("2864165568461618616184867").unwrap().shl(65).to_string(),"105669058452284624487376314655265894269190144");
 
 assert_eq!(Mpz::from_u64(1).shr(0).to_u64().unwrap(), 1);
 assert_eq!(Mpz::from_u64(1).shr(1).to_u64().unwrap(), 0);
 assert_eq!(Mpz::from_u64(1).shr(15).to_u64().unwrap(), 0);
 assert_eq!(Mpz::from_u128(36893488147419103232).shr(65).to_u64().unwrap(),1);
 
 assert_eq!(Mpz::from_u64(16841684685).and(&Mpz::from_u64(1854618418)).to_u64().unwrap(),1786982912);
 assert_eq!(Mpz::from_u64(16841684685).or(&Mpz::from_u64(1854618418)).to_u64().unwrap(),16909320191);
 assert_eq!(Mpz::from_u64(16841684685).xor(&Mpz::from_u64(1854618418)).to_u64().unwrap(),15122337279);
 }
 
