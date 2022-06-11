use std::cmp::Ordering;
use crate::arithmetic::sliceops::*;
use crate::arithmetic::inlineops::*;
use crate::arithmetic::muldiv::*;
use crate::arithmetic::sign::Sign;
use crate::arithmetic::conversion::*;

use crate::primes::PRIMELIST;
use crate::primes::MERSENNE_LIST;

use crate::traits::NumberTheory;



#[derive(Debug,Default,Clone, PartialEq)]
 pub struct Mpz{
      pub(crate)   sign: Sign,
      pub(crate)  limbs: Vec<u64>,
 }
     
 
 impl Mpz{
 
 
   /**
   
   ```
    use number_theory::Mpz;  // includes arbitrary-precision arithmetic
    use number_theory::Sign; // allows sign
    use number_theory::NumberTheory; // includes methods from NumberTheory trait
    
    let bignum = Mpz::new(Sign::Positive, vec![5,5]);
    
    let fortyfour = Mpz::from_u128(44u128);
    
    let fortyfour_neg = Mpz::from_i128(-44i128);
    
    let twopow = Mpz::from_u64(128);
    
    
    assert_eq!("92233720368547758085", bignum.to_string());
    assert_eq!("44", fortyfour.to_string());
    assert_eq!("-44", fortyfour_neg.to_string());
    assert_eq!("128", twopow.to_string());
    
   ```
   
   */
 
  pub fn u_new(limbs: Vec<u64>) -> Self {
         Mpz::from_slice(Sign::Positive, &limbs[..])
  }
  
  pub fn new(sign: Sign, limbs: Vec<u64>)->Self{// New 
         Mpz::from_slice(sign, &limbs[..])
     }
     
  pub fn unchecked_u_new(limbs: Vec<u64>) -> Self{
         Mpz{sign: Sign::Positive, limbs}
  }  
  
  pub fn unchecked_new(sign: Sign, limbs: Vec<u64>) -> Self{
         Mpz{sign, limbs}
  }
  
  pub fn from_slice(sign: Sign, x :&[u64]) -> Self{
             let mut limbs = vec![]; 
             limbs.extend_from_slice(&x[..sig_pos(x)]);
             if limbs.len() == 0{
               limbs.push(0u64)
             }
             Mpz{sign, limbs}
  }   
     
  pub  fn from_u128(x: u128) -> Self{
       let (x_lo,x_hi) = split(x);
       if x_hi == 0 {
        return Mpz::unchecked_new(Sign::Positive, vec![x_lo])
       }
       Mpz::unchecked_new(Sign::Positive, vec![x_lo,x_hi])
 }
 
 pub fn from_i128(x: i128) -> Self{
     if x < 0i128 {
      let (x_lo, x_hi) = split(x.abs() as u128);
      if x_hi == 0 {
      return  Mpz::unchecked_new(Sign::Negative, vec![x_lo])
      }
     return  Mpz::unchecked_new(Sign::Negative, vec![x_lo,x_hi])  
     }
     else{
      return  Mpz::from_u128(x as u128)
     }
 }
 
 pub fn from_u64(x: u64) -> Self{
    Mpz::unchecked_new(Sign::Positive, vec![x])
 }
 
 pub fn from_i64(x: i64) -> Self{
    if x < 0i64 {
     return Mpz::unchecked_new(Sign::Negative, vec![x.abs() as u64])
    }
    return Mpz::unchecked_new(Sign::Positive, vec![x as u64])
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
 
 pub fn rand(len: usize,gen: fn()->u64) -> Self {
      let interim = (0..len).map(|_| rng_64()).collect::<Vec<u64>>();
      Mpz::unchecked_new(Sign::Positive, interim)
 }
 
 pub fn to_string(&self) -> String{
     to_string(self.sign.clone(),self.limbs.clone())
 
 }
    // temporary placeholder function to flip negative zeros
 pub (crate) fn fix_zero(&mut self){
    if self.len() == 0 { self.limbs.push(0u64)}
    if self.len() == 1 && self.limbs[0]== 0 && self.sign == Sign::Negative{
    self.sign = Sign::Positive;
    }
 }
 /**
   Returns the polynomial representation of self in the form of the coefficient vector. Here we see that 50620 to radix-127 is 3x^2 + 17x + 74. 
   Accepts all radix in the interval 0;2^64-1. 
 
 ```
 use number_theory::Mpz;

  let value = Mpz::from_u64(50620);
  assert_eq!(value.to_radix_vec(127),vec![74,17,3])
 ```
 */
 
 pub fn to_radix_vec(&self, radix: u64) -> Vec<u64> {
       let mut k = vec![];
       let mut x = self.clone().limbs;
   loop {
      
      let idx = sig_pos(&x[..]);
      
      if idx == 0 {
        break;
      }
      
      k.push(div_slice(&mut x[..idx],radix));
      
      }
      k.push(x[0]%radix);
      remove_lead_zeros(&mut k);
      return k
 }   
  #[deprecated( note = "u_from_string was originally a quick implementation for usage in other libraries, from_string is stable now so use it instead ~ J.A Sory")]
 pub fn u_from_string(x: &str) -> Option<Self> {
    match from_string(x) {
       Some(y) => Some(Mpz::unchecked_new(Sign::Positive, y)),
       None    => None,
     }
 }
 
 /**
   Conversion from radix-10^n string.  
 
 ```
 use number_theory::Mpz;
  let num = "-3141592653589793238462643383279502884197169399375105820974944592307816406286".to_string();
  let machine = Mpz::from_string(&num).unwrap();
  
  assert_eq!(machine.to_string(),num)
 ```
 */
  
 pub fn from_string(x: &str) -> Option<Self> {
     let ch = x.chars().nth(0).unwrap();
     let mut sign = Sign::Positive;
     let mut k = x;
     if ch == '-'{
       sign = Sign::Negative;
      let mut chars =  x.chars();
      chars.next();
       k = chars.as_str();
     }
     
     if ch == '+'{
       let mut chars = x.chars();
       chars.next();
       k = chars.as_str();
     }
     
     match from_string(k) {
       Some(y) => Some(Mpz::unchecked_new(sign, y)),
       None    => None,
     }
 }
 
  pub  fn zero() -> Self{
         Mpz::unchecked_new(Sign::Positive, vec![0])
  }
  
  pub fn one() -> Self{
         Mpz::unchecked_new(Sign::Positive, vec![1])
  }
     
  pub  fn neg(&mut self){// negation
        self.sign = self.sign.neg();
     }
     
  pub  fn is_one(&self) -> bool{
       if self.len() == 1 && self.limbs[0] == 1{
          return true
       }
         return false
  }   
  
  pub fn is_positive(&self) -> bool{
     self.sign == Sign::Positive
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
 
 pub fn is_mersenne(&self) -> Option<u64> {
    
    let mut flag = 0u64; 
    let mut start = 1u64; 
    let lead = self.limbs.last().unwrap(); 
    for i in 0..64 {
      if *lead == start{
         flag = i;
         break;
      }
      start= (start<<1) + 1;
    }
    if flag == 0 {return None}
    for i in self.limbs[..self.len()-2].iter(){
      if *i != u64::MAX {
      	return None
      }
    }
    
    return Some(flag + 1 + 64u64*(self.len()-1) as u64)
 
 } 
    
 pub  fn set_bit(&mut self, index: usize){ //flips the bit at the index
  
  self.limbs[index/64usize]|=1<<(index%64)
 
 }
 
 pub fn set_sign(&mut self, sign: Sign) {
      self.sign = sign
 } 
     
  pub  fn len(&self)->usize{
           self.limbs.len()
        }
        
  pub fn lead_digit(&self) -> u64{
         *self.limbs[..].last().unwrap()
         
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
       
 pub fn leading_zeros(&self) -> u64 {
    let mut idx : u64 =0;
    
            for i in self.limbs.iter().rev(){
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
              return self.limbs[self.len()-1-idx as usize].leading_zeros() as u64 + 64u64*idx
            }   
 }  
 
 pub fn bit_length(&self) -> u64 {
 	 64u64*self.len() as u64-self.leading_zeros()
 }    
    
 pub fn normalize(&mut self){
      remove_lead_zeros(&mut self.limbs);
      if self.len() == 0{self.limbs.push(0u64)};
 }
 

 /*
 Equality Operations 
 
 */
     
   pub fn u_cmp(&self, other:&Self)-> Ordering{
            cmp_slice(&self.limbs[..],&other.limbs[..])
   } 
 
 
 
  pub fn congruence_u64(&self, n: u64, c: u64) -> bool{
         let mut  interim = mod_slice(&self.limbs[..],n);
         if self.sign == Sign::Negative {
            interim = n - interim;
         }
        interim == c
  }
  
  pub fn add_modinv(&self, n: &Self) -> Self{// additive modular inverse
      let mut k = n.clone();
      let mut selfie = self.clone();

     if self.u_cmp(&n)== Ordering::Greater{
      selfie = selfie.ref_euclidean(&n).1;
     }
     
      sub_slice(&mut k.limbs,&selfie.limbs);
      k.normalize();
      k.sign = Sign::Positive;
      k
  } 

  /**  
  
     ```
           use number_theory::Mpz; 
       let mut one = Mpz::one();
       
       one.successor();  //Applies the successor function ignoring sign
       
       assert_eq!("2", one.to_string())
       
     ```
     
  
  */
  
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
 
  
  
  /*
    Precursor NT functions, unsigned 
  */
     
 pub fn u_quadratic_residue(&self, n: &Self) -> Self{
        self.ref_product(&self).ref_euclidean(n).1
 }   
 
 pub fn u_mul_mod(&self, other: &Self, n: &Self) -> Self{
      self.ref_product(other).ref_euclidean(n).1
  }
   
 pub fn u_mod_pow(&self,  y: &Self, modulo: &Self )->Self{
 
        if modulo == &Mpz::zero(){
          return Mpz::zero()
        }
        
        let mut z = Mpz::one();
        let mut base = self.clone().ref_euclidean(modulo).1;
        let one = Mpz::one();
        

        let mut pow = y.clone();
  
        if pow == Mpz::zero(){
           return z;
        }
        
        if pow == Mpz::one() {
            return self.clone();
        }
        
        while pow.u_cmp(&one) == Ordering::Greater {
       
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
  
  pub fn pow(&self, x: u64)-> Self{ 

  let mut z = Mpz::one();
  let mut base = self.clone();
  let mut pow = x.clone();
  
  if pow ==0 {
    return z 
  }

 while pow > 1 {
  
   if pow%2 == 0 {
      base = base.ref_product(&base);
      pow>>=1;
   }
  
  else{
  
   z = base.ref_product(&z);
   base = base.ref_product(&base);
   pow=(pow-1)>>1;  
   
 }
 }

  base.ref_product(&z)

}


   
 
  
  
  
  pub(crate) fn word_div(&self,x: u64)->(Self, u64){
       let mut quotient = self.clone();
       let remainder = div_slice(&mut quotient.limbs[..], x);
       (quotient,remainder)
  }
  
   fn delta(&self, other: &Self) -> Self{
     if self.u_cmp(other) == Ordering::Greater {
       let mut k = self.clone();
       sub_slice(&mut k.limbs[..],&other.limbs[..]);
       return k 
     }
     if self.u_cmp(other) == Ordering::Less{
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
     k.ref_euclidean(n).1
  }
  
  fn gcd(&self, other: &Self)->Self{
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
  
  pub(crate) fn rho_mpz(&self) -> Self{
   
     let mut x = Mpz::from_u64(2);
     let mut y = Mpz::from_u64(2);
     let mut d = Mpz::one();
  while d == Mpz::one() {
  x = x.mod_sqr_1(&self);
  y = y.mod_sqr_1(&self).mod_sqr_1(&self).ref_euclidean(&self).1 ;
  d = x.delta(&y).gcd(self);
   }
  
    return d
  
  }
  
  /*
 pub fn sqrt(&self) -> Self{
 
 
     let isqrt = |x: u64| { let mut est = x>>((64-x.leading_zeros())/2);
         if est == 0 {
           return 1
         }
    for i in 0..5{
        est = (est + x/est)>>1;
    }
    est
    
    };
 
    if self.len() == 1{
       return Mpz::from_u64(isqrt(self.lead_digit()))
    }
       let two = Mpz::from_u64(2);
    
       let zeros = self.lead_digit().leading_zeros();

       let lead = isqrt((self.lead_digit()>>zeros) + self.limbs[self.len()-2]>>(64-zeros));
       let len = self.bit_length()>>1;
            
           //let checks = 2usize * (self.len() as f64).log2().ceil() as usize + 10usize;
       let mut est = Mpz::from_u64(lead).shl(len as usize); // 30 -1 -1-1-1
       
       for i in 0..35{              // upperbound of 40
         est = ( est.ref_addition(&self.ref_euclidean(&est).0) ).ref_euclidean(&two).0;
       }
       est
    
 
   
 } 
 */

 
 pub fn sqrt(&self) -> Self{
 
    let mut est = self.shr(((self.bit_length()/2)-1) as usize);
   
   let mut count = 0u64;
    loop {
    count +=1;
    let s = est.clone();
    let t = s.ref_addition(&self.euclidean_div(&s).0);
    est = t.shr(1);
    remove_lead_zeros(&mut est.limbs);
    if est.u_cmp(&s) == Ordering::Greater || est.u_cmp(&s) == Ordering::Equal{
      return s
    }
    }
    
    
    
 }
 
 pub fn nth_root(&self, y: u64) -> Self{
      let shift = ((self.bit_length()/y)-1)*(y-1);
      let mut est = self.shr(shift as usize);
      let scalar = Mpz::from_u64(y-1);
      let ymp = Mpz::from_u64(y);
      let mut count = 0u64;
      loop{
      count+=1;
      let s = est.clone();
      let t = s.ref_product(&scalar).ref_addition(&self.euclidean_div(&s.pow(y-1)).0);
      est = t.euclidean_div(&ymp).0;
       if est.u_cmp(&s) == Ordering::Greater || est.u_cmp(&s) == Ordering::Equal{

      return s
    }
      }
 }
 /*
 pub fn nth_root(&self, y: u64) -> Self{
    let nrt = |x: u64, y: u64| { 
    let mut est = x>> ( (y as u32-1)*(64-x.leading_zeros())/y as u32);
         
    for i in 0..5{
        est = ((y-1)*est + x/est.pow(y as u32-1))/y;
    }
    est
    
};

let lead = self.lead_digit();
       let mut est = Mpz::unchecked_new(Sign::Positive, vec![0u64;(self.len()/y as usize)]);
       let root = Mpz::from_u64(y);
       let root_minus =  Mpz::from_u64(y-1);
       est.limbs.push(nrt(lead,y));
       for i in 0..100{
         est = ( (est.ref_product(&root_minus)).ref_addition(&self.ref_euclidean(&est.pow(y-1)).0) ).ref_euclidean(&root).0;
       }
       est

 }
 */
 // Sloppy approximation to natural logarithm
 pub fn ln(&self) -> f64{
     let mut iter_sqrt = 1u64;
     let mut n = self.sqrt(); 
      while n.len() > 1 {
        n = n.sqrt();
        iter_sqrt+=1;
      }
     
      (n.to_u64().unwrap() as f64).ln()*2f64.powf(iter_sqrt as f64)
 }
 
 pub fn log2(&self) -> f64{
     self.ln()*1.4426950408889634
 }
 
 pub fn log10(&self) -> f64{
     self.ln()*0.43429448190325176
 }
 
 pub fn log(&self, log: f64) -> f64{
    self.ln() * log.ln().recip()
 }
 
 pub fn iter_log(&self, log: f64) -> u8 {
      let mut first_log = self.log(log);
      let mut count = 1u8;
      // 1.444667861009766

      while first_log > 1.0 {
         first_log = first_log.log(log) ;
         count+=1
      }
      return count
 }
 
 pub fn eea(&self, other: &Self) -> (Self,Self,Self){
 
       let mut gcd = self.clone();
       let mut new_r = other.clone();
       let mut bezout_1 = Mpz::one();
       let mut new_s = Mpz::zero();
       let mut bezout_2 = Mpz::zero();
       let mut new_t    = Mpz::one();
       
       while new_r != Mpz::zero(){
         let (quo,rem) = gcd.euclidean_div(&new_r);
         let mut temp = new_r.clone();
         new_r = gcd.ref_subtraction(&quo.ref_product(&temp));    
         gcd = temp.clone();
         
         temp = new_s.clone();
         new_s = bezout_1.ref_subtraction(&quo.ref_product(&temp));
         bezout_1 = temp.clone();
         
         temp = new_t.clone();
         new_t = bezout_2.ref_subtraction(&quo.ref_product(&temp));
         bezout_2 = temp.clone();
       
       }
       (gcd,bezout_1,bezout_2) 
 
 }
 
 
 
 /**
   ```
   
      // see the sirp crate for greater extension of this functionality
       use number_theory::Mpz; 
      let factorial_100 = Mpz::sirp(1,100,1,0); 
      
      let doublefact_100 = Mpz::sirp(1,100,2,0);
      assert_eq!(
      "3424322470251197624824643289520818597\
       5118675053719198827915654463488000000000000",doublefact_100.to_string())
   ```
 */
 
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
     /**
     Conditional Interval Product computes the product of integers satisfying an unary function. In this example the unary function is the primality function.
     
   ```
   
       use number_theory::Mpz; 
       use crate::number_theory::NumberTheory;

      let primorial_100 = Mpz::cip(1,100,u64::is_prime); 
      
      
   ```
 */
  
       // conditional interval product
  pub fn cip(infimum: u64, supremum: u64, cond: fn (&u64) -> bool) -> Self{
      let mut cip = Mpz::one();
        let mut acc = 1u64;  // accumulator for factors 
       for i in infimum..supremum+1{
          if cond(&i){
          
            if i >= 4294967296 {  
         acc = i
        }
       if acc < 4294967296{
          acc*=i
       }
       if acc >= 4294967296{
        let mut carry = 0u64;
           carry = scale_slice(&mut cip.limbs[..],acc);

          if carry > 0 {
       
            cip.limbs.push(carry)
      
          }
         acc= 1u64;
       } // end if 
      } // else
     }
     
     let carry = scale_slice(&mut cip.limbs[..],acc);
     if carry > 0{
        cip.limbs.push(carry)
     }
     cip
          }
          
          
  pub fn pi(&self) -> Self{
  
  // let lo = n / (ln - 1. - invln);
     let ln = self.ln();
     let ln_inv = ln.recip();
     let lo_div = Mpz::from_u128( ((ln - 1. - ln_inv).recip()*1E+20) as u128);
     let lo_prod = self.ref_product(&lo_div) ; 
     let div = ln_inv * (1. + ln_inv * (1. + ln_inv * (2. + ln_inv * 7.59)));
     let scalar = Mpz::from_u128((div * 1E+20) as u128);
     let ten = Mpz::from_u64(10).pow(20);
     let prod = self.ref_product(&scalar);
     let total_prod = prod.ref_product(&lo_prod);
      prod.ref_euclidean(&ten).0
  }        
          
          
          
       }
  
  
  
  
  
  
  
