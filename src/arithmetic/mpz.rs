use std::cmp::Ordering;
use crate::arithmetic::sliceops::*;
use crate::arithmetic::inlineops::*;
use crate::arithmetic::muldiv::*;
use crate::arithmetic::sign::Sign;
use crate::arithmetic::conversion::*;

use crate::primes::PRIMELIST;
use crate::primes::MERSENNE_LIST;

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
 
 pub fn rand(len: usize,gen: fn()->u64) -> Self {
      let interim = (0..len).map(|_| rand()).collect::<Vec<u64>>();
      Mpz::new(Sign::Positive, interim)
 }
 
 pub fn to_string(&self) -> String{
     to_string(self.sign.clone(),self.limbs.clone())
 
 }   
  
 pub fn u_from_string(x: &str) -> Option<Self> {
    match from_string(x) {
       Some(y) => Some(Mpz::new(Sign::Positive, y)),
       None    => None,
     }
 }
  
 pub fn from_string(sign: Sign, x: &str) -> Option<Self> {
     match from_string(x) {
       Some(y) => Some(Mpz::new(sign, y)),
       None    => None,
     }
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
     
  pub   fn len(&self)->usize{
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
     
   pub fn cmp(&self, other:&Self)->Ordering{
            cmp_slice(&self.limbs[..],&other.limbs[..])
   } 
 
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
 
 
 
  pub fn congruence_u64(&self, n: u64, c: u64) -> bool{
        mod_slice(&self.limbs[..],n) == c
  }
  
  pub fn add_modinv(&self, n: &Self) -> Self{// additive modular inverse
      let mut k = n.clone();
      let mut selfie = self.clone();

     if self.cmp(&n)== Ordering::Greater{
      selfie = selfie.euclidean(&n).1;
     }
     
      sub_slice(&mut k.limbs,&selfie.limbs);
      k.normalize();
      k.sign = Sign::Positive;
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
  
  /*
      Arithmetic 
  */
  
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
       other.neg();
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
  
  /*
    Precursor NT functions, unsigned 
  */
     
 pub fn u_quadratic_residue(&self, n: &Self) -> Self{
        self.ref_product(&self).euclidean(n).1
 }   
 
 pub fn u_mul_mod(&self, other: &Self, n: &Self) -> Self{
      self.ref_product(other).euclidean(n).1
  }
   
 pub fn u_mod_pow(&self,  y: &Self, modulo: &Self )->Self{
 
        if modulo == &Mpz::zero(){
          return Mpz::zero()
        }
        
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


    
  pub fn sprp(&self, base: Self)->bool{
      let mut p_minus = self.clone();
      let one = Mpz::one();

      sub_slice(&mut p_minus.limbs[..],&one.limbs[..]); //subtract one this will always succeed

      let zeroes = p_minus.trailing_zeros() as usize;

       let d = p_minus.shr(zeroes);  
      let mut x = base.u_mod_pow(&d, self);
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
     let z = Mpz::rand(self.len(),rand).euclidean(&self).1;
     
     if self.sprp(z)==false{ return false} //
   }
   
   return true
  }
  
    // weighted to maintain at most 2^-64 probability of failure, slower than most implementations for small numbers but faster for larger. Values greater than 2^512 receive only two checks, a strong-base 2 and a random-base check. This is due to the fact that the density of pseudoprimes rapidly declines
  pub fn probable_prime(&self) -> bool{
      const CHECK_LUT : [u8;10] = [12,11,9,6,5,5,4,3,2,1];
      let mut check = 1;
     if self.len() < 8{
         check = CHECK_LUT[self.len()]
     }

   let two = Mpz::from_u64(2);

   if self.is_even(){return false}
   if self.is_fermat(){return false}
   
    match self.is_mersenne() {
      Some(x) => {if MERSENNE_LIST.contains(&(x as u32)){return true}
                   else if x < 57885161 {return false}
                   else{return self.llt()} }
      None    => (),
    }
    
   for i in PRIMELIST[1..].iter(){ // apparently optimal on my machine
     if self.congruence_u64(*i as u64,0){
       return false
     }
   }
   
   if self.sprp(two)==false{return false}
 
  let z = self.sprp_check(check as usize);
     
   return z
  
  }
  
  pub fn llt(&self) -> bool{
  	let mut s = Mpz::from_u64(4);
  	let mut p = 0u64;
  	
  	match self.is_mersenne() {
  	  Some(x) => p = x,
  	  None    => return false
  	}
  	
  	for i in 0..(p-2){
  	 // if i%1000==0{println!("{}",i);}
  	  s = s.ref_product(&s);
  	  s.normalize();
  	  sub_slice(&mut s.limbs[..], &[2]);
  	  
  	  s = s.euclidean(&self).1;
  	}
  	s.normalize();
  	if s == Mpz::zero(){return true}
  	return false
  }
  
   /*
   fn lltmod(mut s: Self,m: &Self, p: usize) -> Self {
    while s.bit_length() > p as u64{
       
       s = s.clone().and(&m).addition(s.clone().shr(p)) ;
       
       }
       if &s == m{
         return Mpz::zero();
       }
       return s
   }
  
  pub fn llt(&self) -> bool{
     let mut s = Mpz::from_u64(4);
     let two = Mpz::from_u64(2);
     let p = 521usize;
     'outer :for i in 0..(p-2){
      if i%10==0 {println!("{}",i)}
      s = s.clone().ref_product(&s);
      s.mut_subtraction(two.clone());
      s= Mpz::lltmod(s,self,p);
       
     }
      if s == Mpz::zero(){
         return true
       }
       return false
  }
  */
  // If self is a sophie prime return the safe 
  pub fn is_sophie(&self) -> Option<Self> {
  	
  	if self.is_prime(){
  	  let mut safe = self.shl(1);
  	  let p= safe.clone();
  	  let two = Mpz::from_u64(2);
  	 safe.successor();
  	 
  	 if two.mod_pow(&p,&safe) == Mpz::one() {
  	    return Some(safe)
  	 }
  	 
  	}
  	return None
  }
  
 
  
  
  
  pub(crate) fn word_div(&self,x: u64)->(Self, u64){
       let mut quotient = self.clone();
       let remainder = div_slice(&mut quotient.limbs[..], x);
       (quotient,remainder)
  }
  
   fn delta(&self, other: &Self) -> Self{
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
  
  
  
 pub fn sqrt(&self) -> Self{
        let isqrt = |x: u64| { let mut est = x>>((64-x.leading_zeros())/2);
         
    for i in 0..5{
        est = (est + x/est)>>1;
    }
    est
    
};
       let lead = self.lead_digit();
       let mut est = Mpz::new(Sign::Positive, vec![0u64;(self.len()/2)]);
       let two = Mpz::new(Sign::Positive, vec![2]);
       est.limbs.push(isqrt(lead));
       for i in 0..100{
         est = ( est.addition(self.euclidean(&est).0) ).euclidean(&two).0;
       }
       est
 } 
 
 pub fn nth_root(&self, y: u64) -> Self{
    let nrt = |x: u64, y: u64| { 
    let mut est = x>> ( (y as u32-1)*(64-x.leading_zeros())/y as u32);
         
    for i in 0..5{
        est = ((y-1)*est + x/est.pow(y as u32-1))/y;
    }
    est
    
};

let lead = self.lead_digit();
       let mut est = Mpz::new(Sign::Positive, vec![0u64;(self.len()/y as usize)]);
       let root = Mpz::from_u64(y);
       let root_minus =  Mpz::from_u64(y-1);
       est.limbs.push(nrt(lead,y));
       for i in 0..100{
         est = ( (est.ref_product(&root_minus)).addition(self.euclidean(&est.pow(y-1)).0) ).euclidean(&root).0;
       }
       est

 }
 
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
 /*
 pub fn pi(&self) -> f64{   // fast pi  parallelize probable prime
    let x = self.ln();
    //(self/x )*(1.0 + 1.0/x  + 2.0/(x.ln()*x.ln()))
 }
 */
 /*
 def mod(n,p):
    """ Returns the value of (s**2 - 2) % (2**p -1)"""
    Mp = (1<<p) - 1
    while n.bit_length() > p: # For Python < 2.7 use len(bin(n)) - 2 > p
        n = (n & Mp) + (n >> p)
    if n == Mp:
        return 0
    else:
        return n

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
  
  }
  
  
