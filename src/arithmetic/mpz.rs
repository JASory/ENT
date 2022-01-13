use std::cmp::Ordering;
use crate::arithmetic::sliceops::*;
use crate::arithmetic::inlineops::*;
use crate::arithmetic::muldiv::*;
use crate::arithmetic::sign::Sign;
use crate::arithmetic::conversion::*;

use crate::primes::PRIMELIST;
use crate::traits::NumberTheory;


#[derive(Debug,Clone, PartialEq)]
 pub struct Mpz{
      pub(crate)   sign: Sign,
      pub(crate)  limbs: Vec<u64>,
 }
     
 
 impl Mpz{
  pub fn new(sign: Sign, limbs: Vec<u64>)->Self{// New 
         Mpz{sign, limbs}
     }
     
  pub fn from_u128(x: u128) -> Self{
       let (x_lo,x_hi) = split(x);
       Mpz::new(Sign::Positive, vec![x_lo,x_hi])
 }
 
 pub fn from_i128(x: i128) -> Self{
     if x < 0i128 {
      let (x_lo, x_hi) = split(x.abs() as u128);
     return  Mpz::new(Sign::Negative, vec![x_lo,x_hi])  
     }
     else{
      return  Mpz::from_u128(x as u128)
     }
 }
 
 pub fn from_u64(x: u64) -> Self{
    Mpz::new(Sign::Positive, vec![x])
 }
 
 pub fn to_u64(&self) -> Option<u64>{
       match self.len(){
        0=> Some(0u64),
        1=> Some(self.limbs[0]),
        _=> None,
       }
 }
 
 pub fn to_u128(&self) -> Option<u128>{
       match self.len(){
        0=> Some(0u128),
        1=> Some(self.limbs[0] as u128),
        2=> Some(fuse(self.limbs[1],self.limbs[0])),
        _=> None,
       }
 }
 
 pub fn to_string(&self) -> String{
     to_string(self.sign.clone(),self.limbs.clone())
 
 }   
  
  
 pub from_string(sign: Sign, x: &str) -> Self {
 
 }
 
  pub fn zero() -> Self{
         Mpz::new(Sign::Positive, vec![0])
  }
  
  pub fn one() -> Self{
         Mpz::new(Sign::Positive, vec![1])
  }
     
  pub   fn neg(&mut self){// negation
        self.sign = self.sign.neg();
     }
     
  pub  fn is_even(&self)->bool{// checks if even 
           self.limbs[0]&1==0
       }
       
  pub  fn is_fermat(&self) -> bool{
  
    if (self.limbs[0] != 1) {return false}
    
    let lead =  self.limbs[..].last().unwrap();

    let mut flag = 0u64;
    for i in 0..64{
      if *lead == 1u64<<i{
        flag = 1;   // set flag 
        break;      // end loop
      }
    }
    
    if flag == 0{return false} // if the flag is not set return false
    
    for i in self.limbs[1..self.len()-1].iter(){
     if *i != 0u64 {return false}
    }
    return true 
  }     
    
 pub  fn set_bit(&mut self, index: usize){ //flips the bit at the index
  
  self.limbs[index/64usize]|=1<<(index%64)
 
 } 
     
  pub   fn len(&self)->usize{
           self.limbs.len()
        }
        
        
 pub   fn trailing_zeros(&self)-> u64{// Trailing zeros
          let mut idx : u64 =0;
    
            for i in self.limbs.iter(){
                if i == &0u64{
                   idx+=1;
                }
                else{
                  break;
                }
            }
        if idx == self.len() as u64{
               return 64u64*idx
            }
            else{
              return self.limbs[idx as usize].trailing_zeros() as u64 + 64u64*idx
            }
       }
    
 pub fn normalize(&mut self){
      remove_lead_zeros(&mut self.limbs);
      if self.len() == 0{self.limbs.push(0u64)};
 }
 

 /*
 Equality Operations 
 
 */
     
   pub fn cmp(&self, other:&Self)->Ordering{
            cmp_slice(&self.limbs[..],&other.limbs[..])
   } 
 
 /*
   Shifting operations
 
 */
 
 pub fn mut_shl(&mut self, shift: usize){
    
    let mut k = self.clone();
    let mut trail_zeroes = vec![0;shift/64usize];
 
    let carry = shl_slice(&mut self.limbs[..],(shift%64usize) as u32);
 
 trail_zeroes.extend_from_slice(&self.limbs[..]);
 
 if carry > 0{
    trail_zeroes.push(carry)
 }
 
 self.limbs = trail_zeroes;
 
 }
 
 
 pub fn mut_shr(&mut self, shift: usize){
 
 let mut carry = 0u64;
 
 let mut vector : Vec<u64> = self.limbs.drain((shift/64usize)..self.limbs.len()).collect();
 let sub_shift = shift%64usize;
 
 for i in vector.iter_mut().rev(){
      carry = carry_shr(carry, *i,sub_shift as u32,i);
    
 }
 
 self.limbs = vector;
 }
 
 pub fn shl(&self, shift: usize)->Mpz{
     let mut k = self.clone();
     k.mut_shl(shift);
     k
 }
 
 pub fn shr(&self, shift: usize)->Mpz{
      let mut k = self.clone();
     k.mut_shr(shift);
     k
 }
 
  pub fn congruence_u64(&self, n: u64, c: u64) -> bool{
        mod_slice(&self.limbs[..],n) == c
  }
  
  pub fn add_modinv(&self, n: &Self) -> Self{// additive modular inverse
      let mut k = n.clone();
      sub_slice(&mut k.limbs,&self.limbs);
      k.normalize();
      k
  } 
  
   pub fn successor(&mut self){
     if self.len()==0{self.limbs.push(1)}
     let mut carry = 1u8;
     for i in self.limbs.iter_mut(){
      carry = adc(carry,*i,0,i);
       if carry == 0{
         break;
       }
     }
     if carry > 0u8{
      self.limbs.push(1u64)
     }
  }
  
  pub fn mut_addition(&mut self, mut other: Self){
     let mut carry = 0u8;
 
     if self.sign == other.sign {
        if &self.limbs.len() < &other.limbs.len(){
 
            self.limbs.extend_from_slice(&other.limbs[self.len()..])
        }
     carry = add_slice(&mut self.limbs[..],&other.limbs[..]);
      if carry == 1u8{
        self.limbs.push(1u64)
      }
     }
     
    else{
       if self.cmp(&other)==Ordering::Less{
             carry = sub_slice(&mut other.limbs[..],&self.limbs[..]);
             *self = other;
        }
     else if self.cmp(&other)==Ordering::Equal{
           self.limbs.truncate(0);
           self.limbs.push(0);
           self.sign = Sign::Positive;
        }
       else { 
           sub_slice(&mut self.limbs[..],&other.limbs[..]);
      
        }
    }
    
 }
 
 pub fn addition(&self, other: Self) -> Self{
     let mut k =self.clone();
         k.mut_addition(other);
         k
 }
 
 pub fn mut_subtraction(&mut self, mut other: Self){
       self.neg();
       self.mut_addition(other)
 }
 
 pub fn subtraction(&self, other: Self) -> Self{
       let mut k = self.clone();
       k.mut_subtraction(other);
       k
 }
 
 
  pub fn ref_product(&self, other: &Self) -> Self{
      let mut t = vec![0u64;self.len()+other.len()+1];
     
     mul_slice(&self.limbs[..], &other.limbs[..],&mut t[..]);
     remove_lead_zeros(&mut t);
     Mpz::new(self.sign.mul(&other.sign),t)
  }
  
 pub fn product(&self, other: Self) -> Self{
      let mut t = vec![0u64;self.len()+other.len()+1];
     
     mul_slice(&self.limbs[..], &other.limbs[..],&mut t[..]);
     remove_lead_zeros(&mut t);
     Mpz::new(self.sign.mul(&other.sign),t)
  }
 
  
  
  
  
  pub fn euclidean(&self, other: &Self)->(Self, Self){

    let mut dividend = self.clone();
    
    if dividend == Mpz::zero() { 
        return (Mpz::zero(), Mpz::zero());
    }

    if other.len() == 1 {
        if other.limbs == [1] {
            return (dividend,Mpz::zero());
        }

        let  rem = div_slice(&mut dividend.limbs, other.limbs[0]);

        remove_lead_zeros(&mut dividend.limbs);
        return (dividend, Mpz::new(Sign::Positive, vec![rem]));
    }
     
     
     if dividend.cmp(&other)== Ordering::Equal{
        return (Mpz::one(), Mpz::zero())
     }
     
     if dividend.cmp(&other)== Ordering::Less{
       return (Mpz::zero(), dividend)
     }
   
    let shift = other.limbs.last().unwrap().leading_zeros() as usize;

    if shift == 0 {

      let (quo , rem) = euclidean_slice(&mut dividend.limbs, &other.limbs[..]);
      (Mpz::new(Sign::Positive, quo), Mpz::new(Sign::Positive, rem) )
      
      } 
      
      else {
        let (q, r) = euclidean_slice(&mut dividend.shl(shift).limbs, &other.shl(shift).limbs[..]);
       

       (Mpz::new(Sign::Positive, q), Mpz::new(Sign::Positive, r).shr(shift) )
    }
  }
  
  
     
 pub fn u_quadratic_residue(&self, n: &Self) -> Self{
        self.ref_product(&self).euclidean(n).1
 }   
 
 pub fn u_mul_mod(&self, other: &Self, n: &Self) -> Self{
      self.ref_product(other).euclidean(n).1
  }
   
 pub fn u_mod_pow(&self,  y: &Self, modulo: &Self )->Self{
        let mut z = Mpz::one();
        let mut base = self.clone().euclidean(modulo).1;
        let one = Mpz::one();
        

        let mut pow = y.clone();
        
        if pow ==  Mpz::one(){
           return z
        }
        
        while pow.cmp(&one) == Ordering::Greater {
       
       if pow.len() ==0{
         break;
       }
       
         if pow.is_even(){
           base = base.u_quadratic_residue(&modulo);
           pow.mut_shr(1);
           remove_lead_zeros(&mut pow.limbs);
                 
         }
         
         else{
            z = base.u_mul_mod(&z,&modulo);
            remove_lead_zeros(&mut base.limbs);
            base = base.u_quadratic_residue(&modulo);
            
            sub_slice(&mut pow.limbs[..],&one.limbs[..]);
            pow.mut_shr(1);
            remove_lead_zeros(&mut pow.limbs);

         }
        
        }
        base.u_mul_mod(&z,&modulo)
        //base
        
    }
    
  pub fn sprp(&self, base: u64)->bool{
      let mut p_minus = self.clone();
      let one = Mpz::one();

      sub_slice(&mut p_minus.limbs[..],&one.limbs[..]); //subtract one this will always succeed

      let zeroes = p_minus.trailing_zeros() as usize;

       let d = p_minus.shr(zeroes);  
      let mut x = Mpz::new(Sign::Positive, vec![base]).u_mod_pow(&d, self);
      if x == Mpz::one() || x == p_minus {
         return true
      }
      
      for i in 0..zeroes -1{

        x = x.u_quadratic_residue(&self);

        if x == p_minus{
          return true
        }
      }
      return false
  }
  
  pub fn sprp_check(&self, steps: usize) -> bool{
      if self.len() < 2 {return self.to_u64().unwrap().is_prime()}  // if fits into u64, reduce to 64-bit check 
    
      if self.is_fermat(){return false}
      
      for i in 0..steps{
     let z = rand();
     if self.sprp(z)==false{ return false} //
   }
   
   return true
  }
  
 
  
  
  
  pub(crate) fn word_div(&self,x: u64)->(Self, u64){
       let mut quotient = self.clone();
       let remainder = div_slice(&mut quotient.limbs[..], x);
       (quotient,remainder)
  }
  
  pub fn delta(&self, other: &Self) -> Self{
     if self.cmp(other) == Ordering::Greater {
       let mut k = self.clone();
       sub_slice(&mut k.limbs[..],&other.limbs[..]);
       return k 
     }
     if self.cmp(other) == Ordering::Less{
       let mut k = other.clone();
        sub_slice(&mut k.limbs[..],&self.limbs[..]);
        return k
     }
     else {
       return Mpz::zero()
     }
  }
  
  fn mod_sqr_1(&self, n: &Self) -> Self{
     let mut k = self.ref_product(&self);
     k.successor();
     k.euclidean(n).1
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
  
  pub(crate) fn rho_mpz(&self) -> Self{
     let mut x = Mpz::new(Sign::Positive,vec![2]); let mut y = Mpz::new(Sign::Positive,vec![2]); let mut d = Mpz::one();
  while d == Mpz::one() {
  x = x.mod_sqr_1(&self);
  y = y.mod_sqr_1(&self).mod_sqr_1(&self).euclidean(&self).1 ;
  d = x.delta(&y).gcd(self);
   }
   d
  }
  
 pub fn sirp(infimum: u64, supremum: u64, modulo: u64, residue: u64) -> Self{
 
 let mut sirp = Mpz::one();
  let mut acc = 1u64;  // accumulator for factors 

     for i in infimum..supremum+1{ // inclusive range 
     
       if i % modulo == residue { // if i is of the residue class modulo n then multiply   If you set residue as stop mod n then you get the k-factorials
       
        
        if i >= 4294967296 {
         acc = i
        }
       if acc < 4294967296{
          acc*=i
       }
       if acc >= 4294967296{
        let mut carry = 0u64;
           carry = scale_slice(&mut sirp.limbs[..],acc);

          if carry > 0 {
       
            sirp.limbs.push(carry)
      
          }
         acc= 1u64;
       } // end if 
      } // else
     }
     
     let carry = scale_slice(&mut sirp.limbs[..],acc);
     if carry > 0{
        sirp.limbs.push(carry)
     }
     sirp
     } 
  
  }
  
  
