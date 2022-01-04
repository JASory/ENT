use crate::traits::NumberTheory;
use crate::primes::PRIMELIST;

impl NumberTheory for u8{

 fn is_prime(&self)->bool{
   if *self == 1 || *self == 0 {
     return false
   }
   for i in PRIMELIST[..6].iter(){
    if *self == *i as u8 {return true}
    if self%*i as u8 == 0 {return false}
   }
   return true
 }
 
 fn factor(&self)-> Vec<u8>{
     let mut n = self.clone();
       let twofactors = n.trailing_zeros();
       n>>=twofactors; 
       
       let mut factors : Vec<u8> = vec![];
       
       if twofactors > 0{
          factors.push(2);
          factors.push(twofactors as u8);
       }
      
       for i in PRIMELIST[1..6].iter(){ // strips out small primes
          if n%*i as u8==0{
            factors.push(*i as u8);
            let mut count = 0u8;
          while n%*i as u8==0 {
            count +=1;
            n/=*i as u8;
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
       self.factor().iter().step_by(2).product::<u8>()
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
 
 fn euler_totient(&self) -> u8{
       let factors = self.factor();
       let numerator = factors.iter().step_by(2).map(|x| x -1u8).product::<u8>();
       let denominator = factors.iter().step_by(2).product::<u8>();
       (self/denominator)*numerator
 }
 
 fn quadratic_residue(&self, n: &Self) -> Self{
     ((*self as u16 * *self as u16) % *n as u16) as u8
 }
 
 fn mul_mod(&self, other: &Self, n: &Self) -> Self{
    ((*self as u16 * *other as u16) % *n as u16) as u8
 }
 
 fn mod_pow(&self, p: &Self, modulus: &Self)-> u8{  

  let mut z = 1u16;
  let mut base = *self as u16;
  let n = modulus.clone() as u16;
  let mut pow = p.clone();
  if pow ==0 {
    return z as u8
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

  (base*z % n) as u8

}

}



impl NumberTheory for i8{
  
  fn is_prime(&self) -> bool{
    (self.abs() as u8).is_prime()
  }
  
  fn factor(&self) -> Vec<Self>{
     (self.abs() as u8).factor().iter().map(|x| *x as i8).collect::<Vec<i8>>()
  }
  
  fn radical(&self) -> Self{
     (self.abs() as u8).radical() as i8
  }
  
  fn k_free(&self, k: &Self) -> bool{
      (self.abs() as u8).k_free(&(k.abs() as u8))
  }
  
  fn gcd(&self, other: &Self) -> Self{
      (self.abs() as u8).gcd( &(other.abs() as u8)) as i8
  }
  
  fn euler_totient(&self) -> Self{
     (self.abs() as u8).euler_totient() as i8
  }
  
  fn quadratic_residue(&self, n: &Self) -> Self{
     (self.abs() as u8).quadratic_residue(&(n.abs() as u8)) as i8
  }
  
  fn mul_mod(&self, other: &Self, n: &Self) -> Self{
     let mut a = self.clone();
     let mut b = other.clone();
     let mut modulo = n.abs() ;
     
     if a < 0i8{
        a= modulo + a ;
     }
     if b < 0i8{
        b = modulo + b;
     }
     (a as u8).mul_mod(&(b as u8), &(modulo as u8)) as i8
  }
  
  fn mod_pow(&self, pow: &Self, n: &Self) -> Self{
   let mut a = self.clone();
   if a < 0i8{
      a = n.abs() + self
   }
     (a as u8).mod_pow( &(pow.abs() as u8), &(n.abs() as u8)) as i8
  }
}
