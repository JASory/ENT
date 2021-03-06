use crate::traits::NumberTheory;
use crate::primes::PRIMELIST;
use crate::primes::PRIME_INV_64;
use crate::primes::PRIME_INV_128;

use crate::fjprime64::fjprime_64;
use crate::arithmetic::inlineops::*;


impl NumberTheory for u64{

 fn rng() -> Self {rng_64() }
 
 fn euclidean_div(&self, other: &Self) -> (Self,Self) {
   (*self/ *other, *self%*other)
  }
  
 fn is_prime(&self)->bool{
   if *self < u32::MAX as u64 { // tree down to u32 if it fits
     return (*self as u32).is_prime()
   }
   if *self&1 == 0{return false}
  
  if self < &0x5A2553748E42E8{
  for i in PRIME_INV_64[..].iter(){
  
    if (*self * *i ) < *self{
      return false
    }
  }
  }
  
  if self > &0x5A2553748E42E8{
  for i in PRIME_INV_128[..].iter(){
    if (*self as u128 * *i as u128) < *self as u128 {
      return false
    }
  }
   }
  
     fjprime_64(*self)  
   
 }
 
 fn factor(&self)-> Vec<Self>{
     let mut n = self.clone();
       let twofactors = n.trailing_zeros();
       n>>=twofactors; 
       
       let mut factors : Vec<u64> = vec![];
       
       if twofactors > 0{
          factors.push(2);
          factors.push(twofactors as u64);
       }
      
       for i in PRIMELIST[1..].iter(){ // strips out small primes
          if n% *i as u64==0{
            factors.push(*i as u64);
            let mut count = 0u64;
          while n% *i as u64==0 {
            count +=1;
            n/= *i as u64;
          }
          factors.push(count);
          }
       }
       
       if n == 1 {return factors}
       
        if n.is_prime(){
        factors.push(n);
        factors.push(1);
         return factors
     }
       
       
       
    while n != 1{
          let k = rho_64(n);
           factors.push(k);
           let mut count = 0u64;
      while n%k == 0{
             count+=1;
             n/=k;
           }
           factors.push(count);
       }
       factors
  }
 
 fn radical(&self)-> Self{
       self.factor().iter().step_by(2).product::<u64>()
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
       let numerator = factors.iter().step_by(2).map(|x| x -1u64).product::<u64>();
       let denominator = factors.iter().step_by(2).product::<u64>();
       (self/denominator)*numerator
 }
 
 fn quadratic_residue(&self, n: &Self) -> Self{
                if n == &0 {return 0u64}
     ((*self as u128 * *self as u128) % *n as u128) as u64
 }
 
 fn mul_mod(&self, other: &Self, n: &Self) -> Self{
                if n == &0 {return 0u64}
    ((*self as u128 * *other as u128) % *n as u128) as u64
 }
 
 fn mod_pow(&self, p: &Self, modulus: &Self)-> Self{  
                  if modulus == &0 {return 0u64}
  let mut z = 1u128;
  let mut base = *self as u128;
  let n = modulus.clone() as u128;
  let mut pow = p.clone();
  if pow ==0 {
    return z as u64
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

  (base*z % n) as u64

}

 fn legendre(&self, p: &Self) -> i8 {
    let k = self.mod_pow(&((*p-1)>>1), p);
    if k == 1{return 1};
    if k == *p-1 {return -1};
    return 0
 }
 
 fn checked_legendre(&self, p: &Self) -> Option<i8> {
     if p == &2 || p.is_prime() == false {
          return None
        } 
       Some(self.legendre(&p))
 }
 
 
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

fn checked_jacobi(&self, k: &Self) -> Option<i8>{
    if k > &0 && *k % 2 == 1 {
     return Some(self.jacobi(k))
    }
     return None
 }
 
 fn smooth(&self) -> Self{
     let k = self.factor();
     k[k.len()-2]
 }
 
 fn is_smooth(&self, b: &Self) -> bool{
     &self.smooth() <= b
 }

}



impl NumberTheory for i64{
  
  fn rng() -> Self {rng_64() as i64 }
   
  fn euclidean_div(&self, other: &Self) -> (Self,Self) {
    (*self/ *other, *self%*other)
  }
  
  fn is_prime(&self) -> bool{
    (self.abs() as u64).is_prime()
  }
  
  fn factor(&self) -> Vec<Self>{
     (self.abs() as u64).factor().iter().map(|x| *x as i64).collect::<Vec<i64>>()
  }
  
  fn radical(&self) -> Self{
     (self.abs() as u64).radical() as i64
  }
  
  fn k_free(&self, k: &Self) -> bool{
      (self.abs() as u64).k_free(&(k.abs() as u64))
  }
  
  fn euclid_gcd(&self, other: &Self) -> Self{
      (self.abs() as u64).euclid_gcd(&(other.abs() as u64)) as i64
  }
  
  fn euler_totient(&self) -> Self{
     (self.abs() as u64).euler_totient() as i64
  }
  
  fn quadratic_residue(&self, n: &Self) -> Self{
     (self.abs() as u64).quadratic_residue(&(n.abs() as u64)) as i64
  }
  
  fn mul_mod(&self, other: &Self, n: &Self) -> Self{
     let mut a = self.clone();
     let mut b = other.clone();
     let mut modulo = n.abs() ;
     
     if a < 0i64{
        a= modulo + a ;
     }
     if b < 0i64{
        b = modulo + b;
     }
     (a as u64).mul_mod(&(b as u64), &(modulo as u64)) as i64
  }
  
  fn mod_pow(&self, pow: &Self, n: &Self) -> Self{
   let mut a = self.clone();
   if a < 0i64{
      a = n.abs() + self
   }
     (a as u64).mod_pow( &(pow.abs() as u64), &(n.abs() as u64)) as i64
  }
  
  
   fn legendre(&self, p: &Self) -> i8 {
    let k = self.mod_pow(&((p.abs()-1)>>1), &p.abs());
    if k == 1{return 1};
    if k == p.abs()-1 {return -1};
    return 0
 }
 
  fn checked_legendre(&self, p: &Self) -> Option<i8> {
      if p.abs() == 2 || p.is_prime() == false {
        return None
      } 
     Some(self.legendre(&p))
     
 }
 
 
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

fn checked_jacobi(&self, k: &Self) -> Option<i8>{
    if k > &0 && *k % 2 == 1 {
     return Some(self.jacobi(k))
    }
     return None
 }
 
 fn smooth(&self) -> Self{
     let k = self.factor();
     k[k.len()-2]
 }
 
 fn is_smooth(&self, b: &Self) -> bool{
     self.smooth() <= b.abs()
 }
 
}

  
 fn delta_u64(x: u64, y: u64)->u64{
      if x > y {
          x-y
       }
      else {
          y -x
      }
    }
     // mod sqr plus 1 
    fn mod_sqr1_64(x: u64, n: u64)->u64{
    ((x as u128 * x as u128 + 1 )%n as u128) as u64
   }
 
  // 64-bit pollard rho
  fn rho_64(n: u64)->u64{

  let mut x = 2; let mut y = 2; let mut d = 1;
  
  while d == 1 {
  x = mod_sqr1_64(x,n);
  y = mod_sqr1_64(mod_sqr1_64(y,n),n)%n;
  d = delta_u64(x,y).euclid_gcd(&n)
   }
   d
}

fn mod_sqr(x: &mut u64, y: &mut u64, p: u64){
     *x=((*x as u128 * *x as u128) %p as u128) as u64;
     *y=((*y as u128 * *y as u128) %p as u128) as u64;
    //((x*x % p),(y*y%p) )
}

pub(crate) fn sprp_64(p: u64, base: u64)->bool{// checks if base^p = 1 mod p  or base^(d*2^n)= -1 for some n  
     let zeroes = (p-1).trailing_zeros() as u64; // Breaks number down to p= d*2^n -1
     let d = (p-1)/ (1<<zeroes);
     let mut x = base.mod_pow(&d,&p); // base^d mod p
     if x == 1u64 || x==p-1{   // checks if base^p = 1 mod p  or base^(d*2^n)= -1
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
 
