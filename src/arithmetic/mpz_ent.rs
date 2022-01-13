use crate::NumberTheory;

use crate::arithmetic::sign::Sign;
use crate::arithmetic::mpz::Mpz;
use crate::arithmetic::sliceops::*;
use crate::arithmetic::inlineops::*;
use crate::primes::PRIMELIST;

impl NumberTheory for Mpz{
  
  
  fn quadratic_residue(&self, n: &Self) -> Self{
      let mut p = self.clone();
    if p.sign == Sign::Negative {
      p = p.add_modinv(&n);
    }
     p.u_quadratic_residue(n)
  }
  
  fn mul_mod(&self, other: &Self, n: &Self) -> Self{
        let mut p = self.clone();
        let mut q = other.clone();
    if p.sign == Sign::Negative {
      p = p.add_modinv(n);
    }
    
    if q.sign == Sign::Negative{
       q = q.add_modinv(n)
    }
  
      p.ref_product(&q).euclidean(n).1
      
      
  }
  
  fn mod_pow(&self,  y: &Self, n: &Self )->Self{
        let mut p = self.clone();
    if p.sign == Sign::Negative {
      p = p.add_modinv(&n);
    }
     p.u_mod_pow(y,n)
    }
    
    
  
  fn is_prime(&self) -> bool{
 
    if self.len() < 2 {return self.to_u64().unwrap().is_prime()}  // if fits into u64, reduce to 64-bit check 
    if self.is_even(){return false}
       
    if self.is_fermat(){return false}
    
   for i in PRIMELIST[1..380].iter(){ // apparently optimal on my machine
     if self.congruence_u64(*i as u64,0){
       return false
     }
   }
   
   for i in 0..5{
     let z = rand();
     if self.sprp(z)==false{ return false} //
   }
   
   return true
 
 }
 
 fn gcd(&self, other: &Self)->Self{
         let mut a = self.clone();
         let mut b = other.clone();

         let mut t = Mpz::zero();
        while b != Mpz::zero() {

             t = b.clone();

             b = a.euclidean(&b).1;

             b.normalize();
             a = t;
        }
        return a
    }
    

    
    
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

          
           if n.sprp_check(5){   // stops if prime 
             factors.push(n.clone());
             factors.push(Mpz::one());
             return factors
          }
        //  println!("Next step");
       'outer : while n != Mpz::one(){
          
              let k = n.rho_mpz();
            //  println!("first rho");
              factors.push(k.clone());
              let mut count = 0u64;
             
     'inner : loop {
            let (inner_quo, inner_rem) = n.euclidean(&k);
             
             if inner_rem != Mpz::zero(){
                break 'inner;
             }
             n = inner_quo;
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
   
   (self.euclidean(&denominator).0).ref_product(&numerator)
 
 }
  
  }
