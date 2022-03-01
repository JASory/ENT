use crate::traits::NumberTheory;

use crate::arithmetic::sign::Sign;
use crate::arithmetic::mpz::Mpz;
use crate::arithmetic::sliceops::*;
use crate::arithmetic::inlineops::*;
use crate::primes::PRIMELIST;

use std::cmp::Ordering;

impl NumberTheory for Mpz{
  
  
  fn rng() -> Self {Mpz::unchecked_new(Sign::Positive, vec![rng_64(), rng_64(), rng_64(), rng_64()])}
  
  
  
  
  fn euclidean_div(&self, other:&Self) -> (Self,Self) {
        let (mut quo, mut rem) = self.ref_euclidean(other);
     if self.sign==Sign::Negative && other.sign == Sign::Negative {
        rem.neg();
     }
     else if self.sign==Sign::Positive && other.sign == Sign::Negative {
        quo.neg();
     }
     else if self.sign==Sign::Negative && other.sign == Sign::Positive{
        quo.neg();
        rem.neg();
     }
     
    quo.fix_zero();
    rem.fix_zero();
    (quo,rem)
  }
  
  
  fn quadratic_residue(&self, n: &Self) -> Self{
  
     if n == &Mpz::zero(){
        return Mpz::zero()
     }
     
   let mut p = self.clone();
   
    if p.sign == Sign::Negative {
       p = p.add_modinv(&n);
    }
    
     self.ref_product(&self).ref_euclidean(n).1
  }
  
  
  fn mul_mod(&self, other: &Self, n: &Self) -> Self{
       
       if n == &Mpz::zero(){
        return Mpz::zero()
       }
       
      let mut p = self.clone();
      let mut q = other.clone();
      
    if p.sign == Sign::Negative {
        p = p.add_modinv(n);
    }
    
    if q.sign == Sign::Negative{
       q = q.add_modinv(n)
    }
  
      p.ref_product(&q).ref_euclidean(n).1
    }
  
  fn mod_pow(&self,  y: &Self, n: &Self )->Self{
       if n == &Mpz::zero(){
         return Mpz::zero()
       }
       
     let mut p = self.clone();
     
       if p.sign == Sign::Negative {
          p = p.add_modinv(&n);
       }
       
     p.u_mod_pow(y,n)
    }
    
    
  
  fn is_prime(&self) -> bool{
   if self.probable_prime()==false{return false}
   
   return self.sprp_check(2)    // 2 more strong fermat checks until Lucas test implemented
 
 }
 
 fn euclid_gcd(&self, other: &Self)->Self{
         let mut a = self.clone();
         let mut b = other.clone();

         let mut t = Mpz::zero();
        while b != Mpz::zero() {

             t = b.clone();

             b = a.ref_euclidean(&b).1;

             b.normalize();
             a = t;
        }
        return a
    }
    

    
    //115792089237316195423570985008687907853269984665640564039457584007913129639935
  fn factor(&self) -> Vec<Self>{
     let mut n = self.clone();
     let mut factors : Vec<Self> = vec![];
     let twofactor = n.trailing_zeros();

     if twofactor > 0{
      n.mut_shr(twofactor as usize);
     factors.push(Mpz::from_u64(2u64));
     factors.push(Mpz::from_u64(twofactor));
     }
     
     'outer :  for i in PRIMELIST[1..].iter(){ // skips two as it has already been eliminated
          let (quo,mut rem) = n.word_div(*i as u64);
        
         
          if rem == 0{

           let mut count = 1u64;
             factors.push(Mpz::from_u64(*i as u64));
             n = quo;
            
    'inner : loop {
             
             let (inner_quo, inner_rem) = n.word_div(*i as u64);
             
             if inner_rem != 0{
                break 'inner;
             }
             n = inner_quo;
             count+=1;
            
            }
            
            factors.push(Mpz::from_u64(count)); 
          }

          
          }
          n.normalize();

          if n == Mpz::one() {
            return factors
          }
          
           if n.sprp_check(5){   // stops if prime 
             factors.push(n.clone());
             factors.push(Mpz::one());
             return factors
          }

       'outer : while n != Mpz::one(){
          
              let k = n.rho_mpz();

              factors.push(k.clone());
              let mut count = 0u64;

     'inner : loop {
            let (inner_quo, inner_rem) = n.ref_euclidean(&k);

             if inner_rem != Mpz::zero(){
                break 'inner;
             }
             n = inner_quo;
             n.normalize(); // remove ? 
             count+=1;

           }
           factors.push(Mpz::from_u64(count));

          if n.sprp_check(5){ // stops if  n is prime
            factors.push(n);
            factors.push(Mpz::one());
            break 'outer;
          }
          
          }

          factors
  }
  
  fn radical(&self) -> Self{
  
    let mut rad = Mpz::one();
    
    for i in self.factor().iter().step_by(2){
      rad = rad.ref_product(i)
    }
    
    rad
  }
  
  fn k_free(&self, x: &Self) -> bool{
      
      for i in self.factor()[1..].iter().step_by(2){
        if i == x{
          return false
        }
      }
      return true
  }
  
 fn euler_totient(&self) -> Self{
  
   let mut factors = self.factor();
   
   let mut denominator = Mpz::one();
   let mut numerator = Mpz::one();
   
   for i in factors.iter().step_by(2){
     denominator = denominator.ref_product(i)
   }
   for i in factors.iter_mut().step_by(2){
     sub_slice(&mut i.limbs[..], &[1]);
    numerator = numerator.ref_product(i)
   }
   
   (self.ref_euclidean(&denominator).0).ref_product(&numerator)
 
 }
 
 
 fn legendre(&self, p: &Self) -> i8 {
     let mut p_minus = p.clone();
     p_minus.set_sign(Sign::Positive);
     sub_slice(&mut p_minus.limbs[..], &[1]);
     let pow = p_minus.ref_euclidean(&Mpz::from_u64(2)).0;
     let k = self.mod_pow(&pow,&p);
    if k == Mpz::one() {return 1};
    if k == p_minus {return -1};
    return 0
 }
 
 fn checked_legendre(&self, p: &Self) -> Option<i8> {
     let mut plus = p.clone();
     plus.set_sign(Sign::Positive);
      if plus == Mpz::from_u64(2) {return None}
     match plus.is_prime(){
       true  => Some(self.legendre(&plus)),
       false => None,
     }
 }
 
 
 fn jacobi(&self, k: &Self) -> i8 {
    let mut n = self.clone();
    let mut p = k.clone();
    let mut t = 1i8;
    n = n.euclidean_div(&p).1;

    while n != Mpz::zero(){
     let zeros = n.trailing_zeros(); 
     n.mut_shr(zeros as usize);
     
     if (p.congruence_u64(8,3) || p.congruence_u64(8,5)) && (zeros%2 == 1) { 

            t = -t
     }
    
        std::mem::swap(&mut n, &mut p);
        
        if n.congruence_u64(4,3) && p.congruence_u64(4,3) {

            t = -t;
        }

        n = n.euclidean_div(&p).1;

       
    }
    
    if p == Mpz::one() {
        t
    } 
    
    else {
        0
    }
}

fn checked_jacobi(&self, k: &Self) -> Option<i8>{
    if k.sign == Sign::Positive && k != &Mpz::zero() && k.is_even() == false {
      return Some(self.jacobi(k))
    }
     return None
 }
 
 
  fn smooth(&self) -> Self{
      let k = self.factor();
      k[k.len()-2].clone()
  }

 fn is_smooth(&self,b: &Self) -> bool{
     let mut k = b.clone();
     k.set_sign(Sign::Positive);
     match self.smooth().u_cmp(&k){
      Ordering::Less => true,
      _ => false,
     }
     
     }
  }
