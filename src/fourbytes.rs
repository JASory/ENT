use crate::traits::NumberTheory;
use crate::primes::PRIMELIST;

use crate::fjprime32::fjprime_32;
use crate::arithmetic::inlineops::*;


impl NumberTheory for u32{

 fn rng() -> Self {rng_32() }

 fn euclidean_div(&self, other: &Self) -> (Self,Self) {
   (*self/ *other, *self%*other)
  }

 fn is_prime(&self)->bool{
  
   if *self < u16::MAX as u32 { // tree down to u16 if it fits
     return (*self as u16).is_prime()
   }
   
   for i in PRIMELIST[..30].iter(){
     if *self%*i as u32 == 0{return false}
   }
     fjprime_32(*self)  
   
 }
 
 fn factor(&self)-> Vec<Self>{
     let mut n = self.clone();
       let twofactors = n.trailing_zeros();
       n>>=twofactors; 
       
       let mut factors : Vec<u32> = vec![];
       
       if twofactors > 0{
          factors.push(2);
          factors.push(twofactors as u32);
       }
       
       if n.is_prime(){
        factors.push(n);
        factors.push(1);
         return factors
     }
       
       for i in PRIMELIST[1..].iter(){ // strips out small primes
          if n% *i as u32==0{
            factors.push(*i as u32);
            let mut count = 0u32;
          while n% *i as u32==0 {
            count +=1;
            n/= *i as u32;
          }
          factors.push(count);
          }
       }
       
        if n == 1 {return factors}
       
    while n != 1{
          let k = rho_32(n);
           factors.push(k);
           let mut count = 0u32;
      while n%k == 0{
             count+=1;
             n/=k;
           }
           factors.push(count);
       }
       factors
  }
 
 fn radical(&self)-> Self{
       self.factor().iter().step_by(2).product::<u32>()
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
 
 
 fn euclid_gcd(&self, other: &Self) -> Self{
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
       let numerator = factors.iter().step_by(2).map(|x| x -1u32).product::<u32>();
       let denominator = factors.iter().step_by(2).product::<u32>();
       (self/denominator)*numerator
 }
 
 fn quadratic_residue(&self, n: &Self) -> Self{
           if n == &0 {return 0u32}
     ((*self as u64 * *self as u64) % *n as u64) as u32
 }
 
 fn mul_mod(&self, other: &Self, n: &Self) -> Self{
                if n == &0 {return 0u32}
    ((*self as u64 * *other as u64) % *n as u64) as u32
 }
 
 fn mod_pow(&self, p: &Self, modulus: &Self)-> Self{  
                if modulus == &0 {return 0u32}
  let mut z = 1u64;
  let mut base = *self as u64;
  let n = modulus.clone() as u64;
  let mut pow = p.clone();
  if pow ==0 {
    return z as u32
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

  (base*z % n) as u32

}

 fn legendre(&self, p: &Self) -> i8 {
    let k = self.mod_pow(&((*p-1)>>1), p);
    if k == 1{return 1};
    if k == *p-1 {return -1};
    return 0
 }
 
 fn checked_legendre(&self, p: &Self) -> Option<i8> {
     if *p == 2 {return None}
     match p.is_prime(){
       true  => Some(self.legendre(p)),
       false => None,
     }
 }
 /*
 fn jacobi(&self, k: &Self) -> i8 {
    let mut n = *self;
    let mut p = *k;
    let mut t = 1i8;
    n %= p;
    
    while n != 0 {
     let zeros = n.trailing_zeros(); 
     n>>=zeros;
     
     if (p % 8 == 3 || p % 8 == 5) && (zeros%2 == 1) { 
            t = -t
     }
    
        std::mem::swap(&mut n, &mut p);
        if n % 4 == 3 && p % 4 == 3 {
            t = -t;
        }
        n %= p;
    }
    
    if p == 1 {
        t
    } 
    
    else {
        0
    }
}

fn checked_jacobi(&self, k: &Self) -> i8{
    if k > &0 && *k % 2 == 1 {
       Some(self.jacobi(k))
    }
     return None
 }
 
 */

}



impl NumberTheory for i32{
  
  fn rng() -> Self {rng_32() as i32}
  
  fn euclidean_div(&self, other: &Self) -> (Self,Self) {
   (*self/ *other, *self%*other)
  }
  
  fn is_prime(&self) -> bool{
    (self.abs() as u32).is_prime()
  }
  
  fn factor(&self) -> Vec<Self>{
     (self.abs() as u32).factor().iter().map(|x| *x as i32).collect::<Vec<i32>>()
  }
  
  fn radical(&self) -> Self{
     (self.abs() as u32).radical() as i32
  }
  
  fn k_free(&self, k: &Self) -> bool{
      (self.abs() as u32).k_free(&(k.abs() as u32))
  }
  
  fn euclid_gcd(&self, other: &Self) -> Self{
      (self.abs() as u32).euclid_gcd(&(other.abs() as u32)) as i32
  }
  
  fn euler_totient(&self) -> Self{
     (self.abs() as u32).euler_totient() as i32
  }
  
  fn quadratic_residue(&self, n: &Self) -> Self{
     (self.abs() as u32).quadratic_residue(&(n.abs() as u32)) as i32
  }
  
  fn mul_mod(&self, other: &Self, n: &Self) -> Self{
     let mut a = self.clone();
     let mut b = other.clone();
     let mut modulo = n.abs() ;
     
     if a < 0i32{
        a= modulo + a ;
     }
     if b < 0i32{
        b = modulo + b;
     }
     (a as u32).mul_mod(&(b as u32), &(modulo as u32)) as i32
  }
  
  fn mod_pow(&self, pow: &Self, n: &Self) -> Self{
   let mut a = self.clone();
   if a < 0i32{
      a = n.abs() + self
   }
     (a as u32).mod_pow( &(pow.abs() as u32), &(n.abs() as u32)) as i32
  }
  
  fn legendre(&self, p: &Self) -> i8 {
       //(self.abs() as u32).legendre(&(p.abs() as u32))
    let k = self.mod_pow(&((*p-1)>>1), p);
    if k == 1{return 1};
    if k == *p-1 {return -1};
    return 0
 }
 
  fn checked_legendre(&self, p: &Self) -> Option<i8> {
     (self.abs() as u32).checked_legendre(&(p.abs() as u32))
 }
 
 /*
 fn jacobi(&self, k: &Self) -> i8 {
    let mut n = *self;
    let mut p = *k;
    let mut t = 1i8;
    n %= p;
    
    while n != 0 {
     let zeros = n.trailing_zeros(); 
     n>>=zeros;
     
     if (p % 8 == 3 || p % 8 == 5) && (zeros%2 == 1) { 
            t = -t
     }
    
        std::mem::swap(&mut n, &mut p);
        if n % 4 == 3 && p % 4 == 3 {
            t = -t;
        }
        n %= p;
    }
    
    if p == 1 {
        t
    } 
    
    else {
        0
    }
}

fn checked_jacobi(&self, k: &Self) -> i8{
    if k > &0 && *k % 2 == 1 {
       Some(self.jacobi(k))
    }
     return None
 }
 
 */
}

  // 32-bit pollard rho
 fn delta_u32(x: u32, y: u32)->u32{
      if x > y {
          x-y
       }
      else {
          y -x
      }
    }

  // mod sqr plus 1 
    fn mod_sqr1_32(x: u32, n: u32)->u32{
    ((x as u64 * x as u64 + 1 )%n as u64) as u32
   }
 
  // 64-bit pollard rho
  fn rho_32(n: u32)->u32{

  let mut x = 2; let mut y = 2; let mut d = 1;
  
  while d == 1 {
  x = mod_sqr1_32(x,n);
  y = mod_sqr1_32(mod_sqr1_32(y,n),n)%n;
  d = delta_u32(x,y).euclid_gcd(&n)
   }
   d
}

pub(crate) fn sprp_32(p: u32, base: u32)->bool{// checks if base^p = 1 mod p  or base^(d*2^n)= -1 for some n  
     let zeroes = (p-1).trailing_zeros() as u32; // Breaks number down to p= d*2^n -1
     let d = (p-1)/ (1<<zeroes);
     let mut x = base.mod_pow(&d,&p); // base^d mod p
     if x == 1u32 || x==p-1{   // checks if base^p = 1 mod p  or base^(d*2^n)= -1
       return true
       }
    for _ in 0..zeroes-1{// checks for all d*2^zeroes. One is subtracted since d*2^n was already checked above
     x = x.quadratic_residue(&p);
     if x == p-1 {       // if any d*2^zeroes = p-1  then it passes
       return true
     }
    }
    return false        // otherwise it fails
 }
