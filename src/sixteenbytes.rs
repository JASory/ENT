 
 use crate::primes::PRIMELIST;
 use crate::arithmetic::mpz::Mpz;
 use crate::arithmetic::inlineops::*;
 
 use crate::traits::NumberTheory;
 
 impl NumberTheory for u128 {
 
  fn rng() -> Self {fuse(rng_64(), rng_64()) }
  
  fn euclidean_div(&self, other: &Self) -> (Self,Self) {
   (*self/ *other, *self%*other)
  }
   
  fn is_prime(&self) -> bool{
     if *self < u64::MAX as u128{
       return (*self as u64).is_prime()
     }
     for i in PRIMELIST[..100].iter(){
       if *self%*i as u128 == 0{
         return false
       }
     }
     
     let k = Mpz::from_u128(*self);
     
     return k.sprp_check(5)

  }
  
  fn factor(&self) -> Vec<Self>{
     let mut n = self.clone();
       let twofactors = n.trailing_zeros();
       n>>=twofactors; 
       
       let mut factors : Vec<u128> = vec![];
       
       if twofactors > 0{
          factors.push(2u128);
          factors.push(twofactors as u128);
       }
      
       for i in PRIMELIST[1..].iter(){ // strips out small primes
          if n% *i as u128==0{
            factors.push(*i as u128);
            let mut count = 0u128;
          while n% *i as u128==0 {
            count +=1;
            n/= *i as u128;
          }
          factors.push(count);
          }
       }
       
      if n < u64::MAX as u128{ // 
      
        let large_fact = (n as u64).factor();
        for i in large_fact{
          factors.push(i as u128)
        }
        return factors
      }
       
     else{
        let large_fact = Mpz::from_u128(n).factor();
        
        for i in large_fact{
          factors.push(i.to_u128().unwrap())
        }
        return factors
       
     }  
       
    
  }   
  
  fn radical(&self)-> Self{
       self.factor().iter().step_by(2).product::<u128>()
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
       let numerator = factors.iter().step_by(2).map(|x| x -1u128).product::<u128>();
       let denominator = factors.iter().step_by(2).product::<u128>();
       (self/denominator)*numerator
 }
 
 fn quadratic_residue(&self, n: &Self) -> Self{
               if n == &0 {return 0u128} Mpz::from_u128(*self).u_quadratic_residue(&Mpz::from_u128(*n)).to_u128().unwrap()
 }
 
 fn mul_mod(&self, other: &Self, n: &Self) -> Self{
                if n == &0 {return 0u128}
   Mpz::from_u128(*self).u_mul_mod(&Mpz::from_u128(*other), &Mpz::from_u128(*n)).to_u128().unwrap()
    //((*self as u128 * *other as u128) % *n as u128) as u64
 }
 
 fn mod_pow(&self, p: &Self, modulus: &Self)-> Self{  
                   if modulus == &0 {return 0u128} Mpz::from_u128(*self).u_mod_pow(&Mpz::from_u128(*p),&Mpz::from_u128(*modulus) ).to_u128().unwrap()
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
     return  Some(self.jacobi(k))
    }
     return None
 }
 
 
}  

 impl NumberTheory for i128{
 
 fn rng() -> Self {fuse(rng_64(), rng_64()) as i128}
 
 fn euclidean_div(&self, other: &Self) -> (Self,Self) {
   (*self/ *other, *self%*other)
  }
 
 fn is_prime(&self) -> bool{
    (self.abs() as u128).is_prime()
  }
  
  fn factor(&self) -> Vec<Self>{
     (self.abs() as u128).factor().iter().map(|x| *x as i128).collect::<Vec<i128>>()
  }
  
  fn radical(&self) -> Self{
     (self.abs() as u128).radical() as i128
  }
  
  fn k_free(&self, k: &Self) -> bool{
      (self.abs() as u128).k_free(&(k.abs() as u128))
  }
  
  fn euclid_gcd(&self, other: &Self) -> Self{
      (self.abs() as u128).euclid_gcd(&(other.abs() as u128)) as i128
  }
  
  fn euler_totient(&self) -> Self{
     (self.abs() as u128).euler_totient() as i128
  }
  
  fn quadratic_residue(&self, n: &Self) -> Self{
     (self.abs() as u128).quadratic_residue(&(n.abs() as u128)) as i128
  }
  
  fn mul_mod(&self, other: &Self, n: &Self) -> Self{
     let mut a = self.clone();
     let mut b = other.clone();
     let mut modulo = n.abs() ;
     
     if a < 0i128{
        a= modulo + a ;
     }
     if b < 0i128{
        b = modulo + b;
     }
     (a as u128).mul_mod(&(b as u128), &(modulo as u128)) as i128
  }
  
  fn mod_pow(&self, pow: &Self, n: &Self) -> Self{
   let mut a = self.clone();
   if a < 0i128{
      a = n.abs() + self
   }
     (a as u128).mod_pow( &(pow.abs() as u128), &(n.abs() as u128)) as i128
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
 
 
 
 }
  
  
