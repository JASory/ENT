use crate::{ntrait::NumberTheory,structs::{Certificate,Factorization},result::NTResult};

// Handle u8,u16 and usize values
macro_rules! promoter(
	($($t:ty;$s:ty),* $(,)*) => {$(

   impl NumberTheory for $t{
      
      fn is_unit(&self) -> bool{
          if *self == 1{
            return true
          }
          false
      }

     fn rng() -> $t{
        <$s>::rng() as $t
     }

     fn residue(&self,ring: Self) -> Self{
	(*self as $s).residue(ring as $s) as $t        
     }
     
     fn mul_inverse(&self, ring: Self) -> NTResult<Self>{
        (*self as $s).mul_inverse(ring as $s).map(|n| n as $t)
     }

     fn euclidean_div(&self, other: Self) -> (Self,Self){
        let (quo,rem) = (*self as $s).euclidean_div(other as $s);
        (quo as $t,rem as $t)
     }

     fn fermat(&self, base: Self) -> bool{
        (*self as $s).fermat(base as $s)
     }

     fn strong_fermat(&self, base: Self) -> bool{
	(*self as $s).strong_fermat(base as $s)
     }

     fn is_prime(&self) -> bool{
        (*self as $s).is_prime()
     }	

     fn prime_proof(&self) -> Certificate<Self>{
        let tmp = (*self as $s).prime_proof();
        Certificate::new(tmp.n as $t,tmp.witness as $t,tmp.fctr.iter().map(|x| *x as $t).collect()) 
     }

    fn prime_list(&self,sup: Self) -> Vec<Self>{
       (*self as $s).prime_list(sup as $s).iter().map(|x| *x as $t).collect()
    }

    fn nth_prime(&self) -> NTResult<Self>{
       let res = (*self as $s).nth_prime();
         match res{
            NTResult::Eval(n) => {
            if n > <$t>::MAX as $s{
             return NTResult::Overflow;
            }
            else{
              return NTResult::Eval(n as $t);
            }
            }
            _=> NTResult::Overflow,
          }
    }

    fn pi(&self) -> Self{
        (*self as $s).pi() as $t
    }

    fn prime_gen(x: u32) -> NTResult<Self>{
        if x > Self::BITS-1{
           return NTResult::Overflow;
        }
       <$s>::prime_gen(x).map(|p| p as $t)
    }

    fn factor(&self) -> NTResult<Factorization<Self>>{
       (*self as $s).factor().map(|fctrs| {
         let basevec = fctrs.factor_iter().map(|x| *x as $t).collect();
         return Factorization::from_components(basevec,fctrs.power);
        }
      )
    }

    fn sqrt(&self) -> (Self,Self){
       let (p,q) = (*self as $s).sqrt();
        (p as $t, q as $t)
    }

    fn nth_root(&self, n: Self) -> (Self,Self){
        let (x,y) = (*self as $s).nth_root(n as $s);
        (x as $t, y as $t)
    }

    fn max_exp(&self) -> (Self,Self){
       let (x,y) = (*self as $s).max_exp();
       (x as $t,y as $t)
    }

    fn radical(&self) -> NTResult<Self>{
        let res = (*self as $s).radical();
         match res{
            NTResult::Eval(n) => {
            if n > <$t>::MAX as $s{
             return NTResult::Overflow;
            }
            else{
              return NTResult::Eval(n as $t);
            }
            }
            _=> NTResult::Overflow,
          }
    }


   fn k_free(&self, k: Self) -> bool{
       (*self as $s).k_free(k as $s)
   }

   fn gcd(&self, other: Self) -> Self{
      (*self as $s).gcd(other as $s) as $t
   }
   
   //fn coprime(&self, other: Self) -> bool{
   //    (*self as $s).coprime(other as $s)
  // }

   fn extended_gcd(&self, other: Self) -> (Self,Self,Self){
     let (x,y,g) = (*self as $s).extended_gcd(other as $s);
     (x as $t,y as $t, g as $t)
   }

   fn lcm(&self, other: Self) -> NTResult<Self>{
      let res = (*self as $s).lcm(other as $s);
        match res{
            NTResult::Eval(n) => {
            if n > <$t>::MAX as $s{
             return NTResult::Overflow;
            }
            else{
              return NTResult::Eval(n as $t);
            }
            }
            _=> NTResult::Overflow,
          }
   }
   
   fn euler_totient(&self) -> Self{
       (*self as $s).euler_totient() as $t
   }

   fn jordan_totient(&self, k: Self) -> NTResult<Self>{
       let res = (*self as $s).jordan_totient(k as $s);
        match res{
            NTResult::Eval(n) => {
            if n > <$t>::MAX as $s{
             return NTResult::Overflow;
            }
            else{
              return NTResult::Eval(n as $t);
            }
            }
            _=> NTResult::Overflow,
          }
   }

   fn exponent(&self) -> NTResult<Self>{
      let res = (*self as $s).exponent();
        match res{
            NTResult::Eval(n) => {
            if n > <$t>::MAX as $s{
             return NTResult::Overflow;
            }
            else{
              return NTResult::Eval(n as $t);
            }
            }
            _=> NTResult::Overflow,
          }
   }

   fn dedekind_psi(&self, k: Self) -> NTResult<Self>{
      let res = (*self as $s).dedekind_psi(k as $s);
        match res{
            NTResult::Eval(n) => {
            if n > <$t>::MAX as $s{
             return NTResult::Overflow;
            }
            else{
              return NTResult::Eval(n as $t);
            }
            }
            _=> NTResult::Overflow,
          }
   }

   fn quadratic_residue(&self, n: Self) -> Self{
      (*self as $s).quadratic_residue(n as $s) as $t
   }

   fn checked_quadratic_residue(&self, ring: Self) -> NTResult<Self>{
      let res = (*self as $s).checked_quadratic_residue(ring as $s);
          match res{
            NTResult::Eval(n) => {
            if n > <$t>::MAX as $s{
             return NTResult::Overflow;
            }
            else{
              return NTResult::Eval(n as $t);
            }
            }
            _=> NTResult::Overflow,
          }
   }

   fn product_residue(&self, other: Self, n: Self) -> Self{
      (*self as $s).product_residue(other as $s,n as $s) as $t
   }

   fn checked_product_residue(&self, other: Self, ring: Self) -> NTResult<Self> { 
	let res = (*self as $s).checked_product_residue(other as $s,ring as $s);
	
	    match res{
            NTResult::Eval(n) => {
            if n > <$t>::MAX as $s{
             return NTResult::Overflow;
            }
            else{
              return NTResult::Eval(n as $t);
            }
            }
            _=> NTResult::Overflow,
          }
   }

   fn exp_residue(&self, pow: Self, ring: Self) -> Self{
      (*self as $s).exp_residue(pow as $s,ring as $s) as $t
   }

   fn checked_exp_residue(&self, pow: Self, ring: Self) -> NTResult<Self>{
     let res = (*self as $s).checked_exp_residue(pow as $s, ring as $s);
        match res{
            NTResult::Eval(n) => {
            if n > <$t>::MAX as $s{
             return NTResult::Overflow;
            }
            else{
              return NTResult::Eval(n as $t);
            }
            }
            NTResult::DNE => NTResult::DNE, 
            _=> NTResult::Overflow,
          }
   }

   fn legendre(&self, p: Self) -> i8{
       (*self as $s).legendre(p as $s)
   }

   fn checked_legendre(&self, p: Self) -> NTResult<i8>{
       (*self as $s).checked_legendre(p as $s)
   }

   fn liouville(&self) -> i8{
      (*self as $s).liouville()
   }

   fn derivative(&self) -> NTResult<Self>{
       let res = (*self as $s).derivative();
       match res{
            NTResult::Eval(n) => {
            if n > <$t>::MAX as $s{
             return NTResult::Overflow;
            }
            else{
              return NTResult::Eval(n as $t);
            }
            }
            _=> NTResult::Overflow,
          }
   }

   fn mangoldt(&self) -> f64{
      (*self as $s).mangoldt()
   }

   fn mobius(&self) -> i8{
      (*self as $s).mobius()
   }

   fn jacobi(&self, k: Self) -> i8{
      (*self as $s).jacobi(k as $s)
   }
 
   fn checked_jacobi(&self, k: Self) -> NTResult<i8>{
      (*self as $s).checked_jacobi(k as $s)
   }

   fn kronecker(&self, k: Self) -> i8{
      (*self as $s).kronecker(k as $s)
   }
 
   fn smooth(&self) -> NTResult<Self>{
      let res = (*self as $s).smooth(); 
      match res{
            NTResult::Eval(n) => {
            if n > <$t>::MAX as $s{
             return NTResult::Overflow;
            }
            else{
              return NTResult::Eval(n as $t);
            }
            }
            _=> NTResult::Overflow,
          }
   }

   fn is_smooth(&self, k: Self) -> bool{
     (*self as $s).is_smooth(k as $s)
   }

   fn ord(&self, ring: Self) -> NTResult<Self>{
       let res = (*self as $s).ord(ring as $s); 
       match res{
            NTResult::Eval(n) => {
            if n > <$t>::MAX as $s{
             return NTResult::Overflow;
            }
            else{
              return NTResult::Eval(n as $t);
            }
            }
            _=> NTResult::Overflow,
          }
   }

}
    )*}
);

promoter!(u8;u32,u16;u32,usize;u64);
