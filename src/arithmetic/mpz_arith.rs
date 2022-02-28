use crate::Mpz;
use crate::arithmetic::sliceops::*;
use crate::arithmetic::muldiv::*;
use crate::arithmetic::inlineops::*;
use crate::Sign;
use std::cmp::Ordering;


impl Mpz {


 /*
   Shifting operations
 
 */
 
 
 pub fn mut_and(&mut self, other: &Self){
        for (i,j) in self.limbs.iter_mut().zip(other.limbs.iter()){
          *i = *i&j
        }
 }
 
 pub fn mut_or(&mut self, other: &Self){
        for (i,j) in self.limbs.iter_mut().zip(other.limbs.iter()){
          *i = *i|j
        }
 }
 
 pub fn mut_xor(&mut self, other: &Self){
        for (i,j) in self.limbs.iter_mut().zip(other.limbs.iter()){
          *i = *i^j
        }
 }
 
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
 
 pub fn and(&self, other: &Self) -> Self {
 	let mut k = self.clone();
 	k.mut_and(&other);
 	k
 }
 
 pub fn or(&self, other: &Self) -> Self {
 	let mut k = self.clone();
 	k.mut_or(&other);
 	k
 }
 
 pub fn xor(&self, other: &Self) -> Self {
 	let mut k = self.clone();
 	k.mut_xor(&other);
 	k
 }
 
 
 /*
     Arithmetic operations
 */

 pub fn mut_addition(&mut self, mut other: Self){
     let mut carry = 0u8;
 
     if self.sign == other.sign {
        if &self.limbs.len() < &other.limbs.len(){
            let len = self.len();
            self.limbs.extend_from_slice(&other.limbs[len..]);
            carry = add_slice(&mut self.limbs[..],&other.limbs[..len]);
         
        }
        else{
     carry = add_slice(&mut self.limbs[..],&other.limbs[..]);
      }
      if carry == 1u8{
        self.limbs.push(1u64)
      }
     }
     
    else{
       if self.u_cmp(&other)==Ordering::Less{
             carry = sub_slice(&mut other.limbs[..],&self.limbs[..]);
             *self = other;
             self.normalize();
        }
     else if self.u_cmp(&other)==Ordering::Equal{
           self.limbs.truncate(0);
           self.limbs.push(0);
           self.sign = Sign::Positive;
        }
       else { 
           sub_slice(&mut self.limbs[..],&other.limbs[..]);
           self.normalize();
        }
    }
    
 }
 

 
 pub fn ref_addition(&self, other: &Self) -> Self{
     let mut k =self.clone();
         k.mut_addition(other.clone());
         k
 }
 
 pub fn mut_subtraction(&mut self, mut other: Self){
       other.neg();
       self.mut_addition(other)
 }
 
 pub fn ref_subtraction(&self, other: &Self) -> Self{
       let mut k = self.clone();
       k.mut_subtraction(other.clone());
       k
 }
 

  pub fn ref_product(&self, other: &Self) -> Self{
  
      if self == &Mpz::zero() {
         return Mpz::zero()
      }
      
      if other == &Mpz::zero() {
         return Mpz::zero()
      }
      
      if self.is_one(){
         return Mpz::unchecked_new(self.sign.mul(&other.sign),other.limbs.clone())
      }
      
      if other.is_one(){
         return Mpz::unchecked_new(self.sign.mul(&other.sign),self.limbs.clone())
      }
  
      let mut t = vec![0u64;self.len()+other.len()+1];
     
     mul_slice(&self.limbs[..], &other.limbs[..],&mut t[..]);
     remove_lead_zeros(&mut t);
     Mpz::unchecked_new(self.sign.mul(&other.sign),t)
  }
  
 pub fn mut_product(&mut self, other: Self){
 
      if self == &Mpz::zero() {
         ()
      }
      
      if other == Mpz::zero() {
         self.limbs.truncate(0);
         self.limbs.push(0u64);
         self.sign = Sign::Positive;
      }
      
      if self.is_one(){

         self.limbs.truncate(0);
         self.limbs.extend_from_slice(&other.limbs[..]);
         self.sign = self.sign.mul(&other.sign);
      }
      
      if other.is_one(){
         self.sign = self.sign.mul(&other.sign);
      }
      else {
      let mut t = vec![0u64;self.len()+other.len()+1];
     
     mul_slice(&self.limbs[..], &other.limbs[..],&mut t[..]);
     remove_lead_zeros(&mut t);
     self.limbs = t;
     self.sign = self.sign.mul(&other.sign);
     }
  }
 
  
  
  
  
  pub fn ref_euclidean(&self, other: &Self)->(Self, Self){

    let mut dividend = Mpz::from_slice(Sign::Positive, &self.limbs[..]);
    
    if dividend == Mpz::zero() { 
        return (Mpz::zero(), Mpz::zero());
    }
    
    if dividend.len() == 0usize{
       return (Mpz::zero(),Mpz::zero())
    }
    
    if other.len() == 0usize || other == &Mpz::zero(){
       panic!("Division by zero is undefined in Z")
    }

    if other.len() == 1 {
        if other.limbs == [1] {
            return (dividend,Mpz::zero());
        }

        let  rem = div_slice(&mut dividend.limbs, other.limbs[0]);

        remove_lead_zeros(&mut dividend.limbs);
        return (dividend, Mpz::unchecked_new(Sign::Positive, vec![rem]));
    }
     
     
     if dividend.u_cmp(&other)== Ordering::Equal{
        return (Mpz::one(), Mpz::zero())
     }
     
     if dividend.u_cmp(&other)== Ordering::Less{
       return (Mpz::zero(), dividend)
     }
   
    let shift = other.limbs.last().unwrap().leading_zeros() as usize;

    if shift == 0 {

      let (quo , rem) = euclidean_slice(&mut dividend.limbs, &other.limbs[..]);
      (Mpz::unchecked_new(Sign::Positive, quo), Mpz::unchecked_new(Sign::Positive, rem) )
      
      } 
      
      else {
        let (q, r) = euclidean_slice(&mut dividend.shl(shift).limbs, &other.shl(shift).limbs[..]);
       

       (Mpz::unchecked_new(Sign::Positive, q), Mpz::unchecked_new(Sign::Positive, r).shr(shift) )
    }
  }
}
