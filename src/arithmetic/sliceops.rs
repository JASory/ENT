use crate::arithmetic::inlineops::*;
use crate::arithmetic::sign::Sign;
use std::cmp::Ordering;

/*
    Slice operations  
*/

   // Compare slices 
pub(crate) fn cmp_slice(x: &[u64], y: &[u64])->Ordering{
 
      if x.len() > y.len(){
          return Ordering::Greater
          }
      if x.len() < y.len(){
          return Ordering::Less
          }
      else{
          Iterator::cmp(x.iter().rev(), y.iter().rev())
          }
  }

/*
   Bitwise slice operations 
*/

   // slice shift-right
pub(crate) fn shr_slice(x: &mut[u64],shift: u32)-> u64{

  let mut carry= 0u64;
  
  for i in x.iter_mut().rev(){
  
   carry = carry_shr(carry, *i, shift,i);
  
  }
  
  carry

}  
  // slice shift-left
pub(crate) fn shl_slice(x: &mut[u64], shift: u32)-> u64{

   let mut carry = 0u64;

     for i in x.iter_mut(){
         
           carry = carry_shl(carry, *i, shift,i);
           
      }
      
   carry
}

/*

 Arithmetic slice operations 
 
*/
#[inline]
pub(crate) fn add_slice(x: &mut [u64],y: &[u64])->u8{  //first number must be larger

    let mut carry = 0u8;
 
    let (lo,hi) = x.split_at_mut(y.len()); //split to make lo equal in length to y
    
    for (i,j) in lo.iter_mut().zip(y.iter()){                               //add equal 
            carry = adc(carry,*i,*j,i);
        }
      
      if carry > 0u8{// if carry is greater than zero, propagate it through the rest of the array until there is no carry left
       
       for k in hi.iter_mut(){   //add the carry until there is no carry
       carry = adc(carry,*k,0u64,k);
       if carry == 0u8{
       break;
       }
       }
      
      }
   carry
 }
 
 
 /*
    Subtraction x-y
 */
 
pub(crate) fn sub_slice(x: &mut [u64],y: &[u64])->u8{  //first number must be larger

 let mut carry = 0u8;
 
 let (lo,hi) = x.split_at_mut(y.len()); //split to make equal
    
 for (i,j) in lo.iter_mut().zip(y.iter()){                               //add equal 
        carry = sbb(carry,*i,*j,i);
       }
      
      if carry > 0u8{
       
       for k in hi.iter_mut(){   //add the carry until there is no carry
       carry = sbb(carry,*k,0u64,k);
       if carry == 0u8{
       break;
       }
       }
      
      }
   carry

}

pub (crate) fn  scale_slice(x: &mut [u64], scale : u64)->u64{
      
   let mut carry = 0u64;
   
   for i in x.iter_mut() {
     carry = carry_mul(carry,*i,scale, i);
   }   
   carry
  }

  // divides inplace and returns remainder
pub(crate) fn div_slice(x : &mut[u64], divisor: u64 )->u64{
 
 let mut carry = 0u64;
   
   for i in x.iter_mut().rev() {
     carry = carry_div(carry,*i,divisor, i);
   }   
   carry
 }
 
  // divides inplace and returns remainder
pub(crate) fn mod_slice(x : &[u64], divisor: u64 )->u64{
 
 let mut carry = 0u64;
   
   for i in x.iter().rev() {
     carry = carry_mod(carry,*i,divisor);
   }   
   carry
 }
