use crate::traits::NumberTheory;
use crate::primes::PRIMELIST;
/*
pub const PRIMELIST: [u16; 54] = [// list of all primes less than 2^8
           2,  3,  5,   7,  11,  13,  17,  19,  23,  29,  31,
          37, 41, 43,  47,  53,  59,  61,  67,  71,  73,  79, 
          83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 
         139,149,151, 157, 163, 167, 173, 179, 181, 191, 193,  //54
         197,199,211, 223, 227, 229, 233, 239, 241, 251
        ];
*/

impl NumberTheory for u16{

 fn is_prime(&self)->bool{
   if *self == 1 || *self == 0 {
     return false
   }
   if *self < u8::MAX as u16 { // tree down to u8 if it fits
     return (*self as u8).is_prime()
   }
        
   for i in PRIMELIST[..54].iter(){
    if self == i {return true}
    if self%i == 0 {return false}
   }
   return true
 }
 
 fn factor(&self)-> Vec<Self>{
     let mut n = self.clone();
       let twofactors = n.trailing_zeros();
       n>>=twofactors; 
       
       let mut factors : Vec<u16> = vec![];
       
       if twofactors > 0{
          factors.push(2);
          factors.push(twofactors as u16);
       }
       
       for i in PRIMELIST[1..54].iter(){ // strips out small primes
          if n%i==0{
            factors.push(*i);
            let mut count = 0u16;
          while n%i==0 {
            count +=1;
            n/=i;
          }
          factors.push(count);
          }
       }
       if n > 1 {
        factors.push(n);
        factors.push(1);
       }
       factors  
 }
 
 fn radical(&self)-> Self{
       self.factor().iter().step_by(2).product::<u16>()
   }
   
 fn k_free(&self, k: &Self)->bool{
       let factors = self.factor();
      for i in 0..factors.len(){
        if factors[i] == *k{
          if i == 0{
            ();
          }
          else{
            return false
          }
        }
      }
      return true
   }
 
 
 fn gcd(&self, other: &Self) -> Self{
     let mut a = self.clone();
     let mut b = other.clone();
     if b == 0 
    { return a; } 

    else if a == 0
    { return b; }

      let self_two_factor = a.trailing_zeros();  
      let other_two_factor = b.trailing_zeros(); 
      let min_two_factor = std::cmp::min(self_two_factor, other_two_factor);
        a >>= self_two_factor;
	b >>= other_two_factor;
      loop {
         
         if b > a {
             std::mem::swap(&mut b, &mut a);
         }
         a -= b;

         if a == 0 {
             return b << min_two_factor;
         }
         a >>= a.trailing_zeros();
     }
 }
 
 fn euler_totient(&self) -> Self{
       let factors = self.factor();
       let numerator = factors.iter().step_by(2).map(|x| x -1u16).product::<u16>();
       let denominator = factors.iter().step_by(2).product::<u16>();
       (self/denominator)*numerator
 }
 
 fn quadratic_residue(&self, n: &Self) -> Self{
     ((*self as u32 * *self as u32) % *n as u32) as u16
 }
 
 fn mul_mod(&self, other: &Self, n: &Self) -> Self{
    ((*self as u32 * *other as u32) % *n as u32) as u16
 }
 
 fn mod_pow(&self, p: &Self, modulus: &Self)-> Self{  

  let mut z = 1u32;
  let mut base = *self as u32;
  let n = modulus.clone() as u32;
  let mut pow = p.clone();
  if pow ==0 {
    return z as u16
  }

 while pow > 1 {
  
   if pow%2 == 0 {
      base = base*base % n ;
      pow>>=1;
   }
  
  else{
  
   z = base*z % n;
   base = base*base % n;
   pow=(pow-1)>>1;  
   
 }
 }

  (base*z % n) as u16

}

}



impl NumberTheory for i16{
  
  fn is_prime(&self) -> bool{
    (self.abs() as u16).is_prime()
  }
  
  fn factor(&self) -> Vec<Self>{
     (self.abs() as u16).factor().iter().map(|x| *x as i16).collect::<Vec<i16>>()
  }
  
  fn radical(&self) -> Self{
     (self.abs() as u16).radical() as i16
  }
  
  fn k_free(&self, k: &Self) -> bool{
      (self.abs() as u16).k_free(&(k.abs() as u16))
  }
  
  fn gcd(&self, other: &Self) -> Self{
      (self.abs() as u16).gcd( &(other.abs() as u16)) as i16
  }
  
  fn euler_totient(&self) -> Self{
     (self.abs() as u16).euler_totient() as i16
  }
  
  fn quadratic_residue(&self, n: &Self) -> Self{
     (self.abs() as u16).quadratic_residue(&(n.abs() as u16)) as i16
  }
  
  fn mul_mod(&self, other: &Self, n: &Self) -> Self{
     let mut a = self.clone();
     let mut b = other.clone();
     let mut modulo = n.abs() ;
     
     if a < 0i16{
        a= modulo + a ;
     }
     if b < 0i16{
        b = modulo + b;
     }
     (a as u16).mul_mod(&(b as u16), &(modulo as u16)) as i16
  }
  
  fn mod_pow(&self, pow: &Self, n: &Self) -> Self{
   let mut a = self.clone();
   if a < 0i16{
      a = n.abs() + self
   }
     (a as u16).mod_pow( &(pow.abs() as u16), &(n.abs() as u16)) as i16
  }
}
