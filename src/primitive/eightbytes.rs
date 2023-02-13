use crate::traits::NumberTheory;
use crate::traits::Reduction;
use crate::result::NTResult;

use crate::arithmetic::inlineops::*;
use crate::montgomery::*;

use crate::data::primes::PRIMELIST;
use crate::data::primes::PRIME_INV_128;
use crate::data::primes::PRIME_INV_64;

use crate::data::hashtable::BASE_33;
use crate::data::hashtable::BASE_34;
use crate::data::hashtable::BASE_35;

use crate::data::hashtable::BASES_35_64;


fn detect_pseudo(x: u64) -> bool {
    for i in 2..16 {
        let sq = (x - 1) / i;
        let k = (sq as f64).sqrt() as u64;
        if ((k * k + k).wrapping_mul(i)).wrapping_add( k + 1) == x {
            return true;
        }
    }
    return false;
}

impl NumberTheory for u64 {
    fn rng() -> Self {
        rng_64()
    }

    fn residue(&self, ring: &Self) -> Self{
      if ring == &0{
        return *self
      }
        *self % *ring
    }
     
    fn euclidean_div(&self, other: &Self) -> (Self, Self) {
        (*self / *other, *self % *other)
    }

    fn is_sprp(&self, base: &Self) -> bool {
        sprp_64(*self, *base)
    }

    fn is_prime(&self) -> bool {
        if *self < u32::MAX as u64 {
            // tree down to u32 if it fits
            return (*self as u32).is_prime();
        }
        if *self & 1 == 0 {
            return false;
        }

        if self < &0x5A2553748E42E8 {
            for i in PRIME_INV_64[..256].iter() {
                if ((*self).wrapping_mul(*i)) < *self {
                    return false;
                }
            }
        }

        if self > &0x5A2553748E42E8 {
            for i in PRIME_INV_128[..128].iter() {
                if ((*self as u128).wrapping_mul(*i)) < *self as u128 {
                    return false;
                }
            }
        }

        // TODO : Reduce to a single hashtable if possible to minimize if branching

        if *self < 8589934592 {
            return self.is_sprp(&(BASE_33[((*self ^ 0x100000000) >> 24) as usize] as u64));
        }

        if *self < 17179869184 {
            return self.is_sprp(&(BASE_34[((*self ^ 0x200000000) >> 24) as usize] as u64));
        }
        if *self < 34359738368 {
            return self.is_sprp(&(BASE_35[((*self ^ 0x400000000) >> 25) as usize] as u64));
        }

        if !self.is_sprp(&2) {
            return false;
        }

        if detect_pseudo(*self) {
            return false;
        }

        let idx = ((*self as u32).wrapping_mul(3301793688) >> 17) as usize;

        self.is_sprp(&(BASES_35_64[idx] as u64))
    }

    fn prime_proof(&self) -> (bool, Vec<Self>) {
        if *self == 2 {
            return (true, vec![3]);
        }

        let x_minus = *self - 1;
        let fctrs = x_minus
            .factor()
            .iter()
            .step_by(2)
            .map(|y| *y)
            .collect::<Vec<Self>>();

        let mut certificate = vec![2];

        certificate.extend_from_slice(&fctrs[..]);

        loop {
            // loops until it has either been shown to be prime or composite

            let mut witness = Self::rng() % (self - 2) + 2;

            'witness: loop {
                if witness.gcd(&self) == 1 {
                    break 'witness;
                }
                witness += 1;
            }

            if witness.exp_residue(&x_minus, &self) != 1 {
                // If any witness says it's composite then it is
                certificate[0] = witness;

                return (false, certificate);
            }

            'inner: for (idx, j) in fctrs.iter().enumerate() {
                if witness.exp_residue(&((*self - 1) / j), self) == 1 {
                    break 'inner;
                }
                if idx == fctrs.len() - 1 {
                    certificate[0] = witness;

                    return (true, certificate);
                }
            }
        }
    }

    fn prime_list(&self, sup: &Self) -> Vec<Self> {
        let inf = std::cmp::min(*self, *sup);
        let mut hi = std::cmp::max(*self, *sup);

        if hi < u64::MAX {
            hi += 1;
        }

        let mut primevector = vec![];

        for i in inf..hi {
            if i.is_prime() {
                primevector.push(i)
            }
        }
        primevector
    }

    fn nth_prime(&self) -> NTResult<Self> {
        let mut count = 0u64;
        let mut start = 0u64;
        
        if *self < 203280221{
          return (*self as u32).nth_prime().map(|y| y as u64)
        }
        
        if *self > 425656284035217743{
          return NTResult::Overflow
        }
    
        loop {
            start += 1;

            if start == Self::MAX {
                return NTResult::Overflow;
            }
            if start.is_prime() {
                count += 1;
            }
            if count == *self {
                return NTResult::Eval(start);
            }
        }
    }

    fn pi(&self) -> Self {
    if  self.reducible(){
         return (*self as u32).pi() as u64
       }
        let mut count = 0u64;
        for i in 0u64..*self {
            if i.is_prime() {
                count += 1;
            }
        }
        count
    }

    fn prime_gen(k: u32) -> NTResult<Self> {
        if k > 64 {
           return NTResult::Overflow;
        }
        if k < 33 {
            return u32::prime_gen(k).map(|x| x as u64);
        }
        let form = (1 << (k - 1)) + 1;
        let bitlength = form - 2;

        loop {
            let p = u64::rng();
            if ((p & bitlength) | form).is_prime() {
                return NTResult::Eval((p & bitlength) | form);
            }
        }
    }

    fn factor(&self) -> Vec<Self> {
        if self < &4294967295{
          return (*self as u32).factor().iter().map(|x| *x as u64).collect::<Vec<u64>>()
        }
        
        let mut n = *self;
        let twofactors = n.trailing_zeros();
        n >>= twofactors;

        let mut factors: Vec<u64> = vec![];

        if twofactors > 0 {
            factors.push(2);
            factors.push(twofactors as u64);
        }

        for i in PRIMELIST[1..].iter() {
            // strips out small primes
            if n % *i as u64 == 0 {
                factors.push(*i as u64);
                let mut count = 0u64;
                while n % *i as u64 == 0 {
                    count += 1;
                    n /= *i as u64;
                }
                factors.push(count);
            }
        }

        if n == 1 {
            return factors;
        }

        if n.is_prime() {
            factors.push(n);
            factors.push(1);
            return factors;
        }

        while n != 1 {
            let k = rho_64(n);
            factors.push(k);
            let mut count = 0u64;
            while n % k == 0 {
                count += 1;
                n /= k;
            }
            factors.push(count);
        }
        factors
    }
    
    fn checked_factor(&self) -> NTResult<Vec<Self>>{
     
     if *self == 0{
       return NTResult::InfiniteSet
     }
     if *self == 1{
       return NTResult::DNE
     }
     
     NTResult::Eval(self.factor())
    
    }


    fn sqrt(&self) -> (Self, Self) {
        if *self < 0x100000000 {
            return ((*self as u32).sqrt().0 as u64, 0);
        }
        let mut est = (*self as f64).sqrt() as Self + 1;

        loop {
            let s = est;
            let t = s + *self / s;
            est = t >> 1;
            if est >= s {
                return (s, 0);
            }
        }
    }

    fn nth_root(&self, n: &Self) -> (Self, Self) {
        if *n > 63 {
            return (1, 0);
        }

        if *n == 1 {
            return (*n, 0);
        }

        if *n == 0 {
            panic!("No integer is a zeroth factor ")
        }

        let mut est = (*self as f64).powf((*n as f64).recip()) as Self + 1;

        loop {
            let s = est;
            let t = (*n - 1) * s + *self / s.pow(*n as u32 - 1);
            est = t / *n;
            if est >= s {
                return (s, 0);
            }
        }
    }

    fn max_exp(&self) -> (Self,Self){
    
      for i in 1..64{
      let p = 64-i;
      let base = self.nth_root(&p).0;
         if base.pow(p as u32) == *self{
           return(base,p)
         }
      }
     return (*self,1)
    }    
    
    fn radical(&self) -> NTResult<Self> {
      if self.reducible(){
         return (*self as u32).radical().map(|kishum| kishum as u64)
      }
        self.checked_factor().map(|y| y.iter().step_by(2).product::<Self>())
    }

    fn k_free(&self, k: &Self) -> bool {
       if self.reducible(){
         return (*self as u32).k_free(&(*k as u32))
      }
        let factors = self.factor();
        for (idx, el) in factors.iter().enumerate() {
            if el == k && idx != 0 {
                return false;
            }
        }
        true
    }

    fn gcd(&self, other: &Self) -> Self {
        let mut a = *self;
        let mut b = *other;
        if b == 0 {
            return a;
        } else if a == 0 {
            return b;
        }

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

    fn extended_gcd(&self, other: &Self) -> (Self, Self, Self) {
        let mut gcd: u64 = *self;
        let mut new_r: u64 = *other;
        let mut bezout_1: u64 = 1;
        let mut new_s: u64 = 0;
        let mut bezout_2: u64 = 0;
        let mut new_t: u64 = 1;

        while new_r != 0 {
            let quotient = gcd / new_r;
            let mut temp: u64 = new_r;
            new_r = gcd - quotient * temp;
            gcd = temp;

            temp = new_s;
            if bezout_1 < quotient.product_residue(&temp, other) {
                // First bezout coefficient is computed over Z[q]
                new_s = *other - (quotient.product_residue(&temp, other) - bezout_1)
            } else {
                new_s = bezout_1.wrapping_sub(quotient * temp);
            }

            bezout_1 = temp;

            temp = new_t;
            if bezout_2 < quotient.product_residue(&temp, self) {
                // Second bezout coefficient is computed over Z[p]
                new_t = *self - (quotient.product_residue(&temp, self) - bezout_2)
            } else {
                new_t = bezout_2.wrapping_sub(quotient.product_residue(&temp, self));
            }
            bezout_2 = temp
        }
        (gcd, bezout_1, bezout_2)
    }

    fn lcm(&self, other: &Self) -> Self {
       if self == &0 && other == &0{
           return 0
         }
        let cf = self.gcd(other);
        (*self / cf) * *other
    }

    fn checked_lcm(&self, other: &Self) -> NTResult<Self> {
        if self == &0 && other == &0{
           return NTResult::Eval(0)
         }
        let cf = self.gcd(other);
        let (v, flag) = (*self / cf).overflowing_mul(*other);
        if flag {
            return NTResult::Overflow;
        }
        NTResult::Eval(v)
    }

    fn euler_totient(&self) -> Self {
      if *self < u32::MAX as u64{
        return (*self as u32).euler_totient() as u64
      }
        let factors = self.factor();
        let numerator = factors.iter().step_by(2).map(|x| x - 1u64).product::<u64>();
        let denominator = factors.iter().step_by(2).product::<u64>();
        (self / denominator) * numerator
    }

    fn jordan_totient(&self, k: &Self) -> NTResult<Self> {
        if *k > u32::MAX as u64 {
            return NTResult::Overflow;
        }
        
        if *self < 2{
           return NTResult::Eval(*self)
        }
        
        let (coef, flag) = self.overflowing_pow(*k as u32);

        if flag {
            return NTResult::CompOverflow;
        }

        let mut denom = 1u64;
        let mut numer = 1u64;

        for i in self.factor().iter().step_by(2) {
            let pow = i.pow(*k as u32);

            denom = denom * pow;

            numer *= pow - 1;
        }

        NTResult::Eval(numer * (coef / denom))
    }
    
    fn carmichael_totient(&self) -> NTResult<Self>{ 
    
       if  *self < 4294967295{
         return (*self as u32).carmichael_totient().map(|x| x as u64)
       }
       
       let fctr = self.factor();
       let base = fctr.iter().step_by(2).map(|z| *z).collect::<Vec<Self>>();
       let mut result = 1;
      for (idx,el) in base.iter().enumerate(){
        if el == &2 && fctr[1] > 2{
         let phi = ((el.pow(fctr[2*idx+1] as u32) /el) *(el-1)) /2;
          result = result.lcm(&phi);
        }
       else{
         let phi =  (el.pow(fctr[2*idx+1] as u32)/el)*(el-1);
         result = result.lcm(&phi);
       } 
      }
     NTResult::Eval(result)
    }

    fn dedekind_psi(&self, k: &Self) -> NTResult<Self> {
        if *self == 0{
         return NTResult::Infinite
       }
       
        let (k2, flag) = k.overflowing_shl(1);
        if flag {
            return NTResult::Overflow;
        }
        self.jordan_totient(&k2).map(|y| y/self.jordan_totient(k).unwrap())
    }

   fn quadratic_residue(&self, n: &Self) -> Self {
        if n == &0 {
            return self.wrapping_mul(*self)
        }
        ((*self as u128 * *self as u128) % *n as u128) as Self
    }
    
    fn checked_quadratic_residue(&self, n: &Self) -> NTResult<Self> {
        if n == &0 {
            return NTResult::from_option(self.checked_mul(*self),NTResult::Overflow)
        }
        NTResult::Eval(((*self as u128 * *self as u128) % *n as u128) as Self)
    }
    

    fn product_residue(&self, other: &Self, n: &Self) -> Self {
        if n == &0 {
            return self.wrapping_mul(*other)
        }
        ((*self as u128 * *other as u128) % *n as u128) as Self
    }
    
    fn checked_product_residue(&self, other: &Self, n: &Self) -> NTResult<Self> {
        if n == &0 {
            return NTResult::from_option(self.checked_mul(*self),NTResult::Overflow)
        }
        NTResult::Eval(((*self as u128 * *other as u128) % *n as u128) as Self)
    }

    fn exp_residue(&self, p: &Self, modulus: &Self) -> Self {
        if modulus == &0 {
            return self.pow(*p as u32);
        }

        pow_64(*self, *p, *modulus)
    }

    fn checked_exp_residue(&self, p: &Self, modulus: &Self) -> NTResult<Self> {
        if modulus == &0 {
            if *p > u32::MAX as u64 {
                return NTResult::Overflow;
            }
            match self.checked_pow(*p as u32) {
                Some(x) => return NTResult::Eval(x),
                None => return NTResult::Overflow,
            };
        }

        NTResult::Eval(pow_64(*self, *p, *modulus))
    }

    fn legendre(&self, p: &Self) -> i8 {
        let k = self.exp_residue(&((*p - 1) >> 1), p);
        if k == 1 {
            return 1;
        };
        if k == *p - 1 {
            return -1;
        };
        0i8
    }

    fn checked_legendre(&self, p: &Self) -> NTResult<i8> {
        if p == &2 || !p.is_prime() {
            return NTResult::Undefined;
        }
        NTResult::Eval(self.legendre(p))
    }

    fn liouville(&self) -> i8 {
      if self.reducible(){
        return (*self as u32).liouville()
      }
        let primeomega = self.factor()[1..].iter().step_by(2).sum::<Self>();
        if primeomega & 1 == 0 {
            return 1;
        }
        return -1;
    }
    
    fn derivative(&self) -> NTResult<Self> {
       if *self < 94 {
         return (*self as u8).derivative().map(|y| y as Self)
       }
       let fctr = self.factor();
       let mut sum : u64 = 0;
       
     for i in 0..fctr.len() / 2 {
        match sum.checked_add(fctr[2 * i + 1] * (*self / fctr[2 * i])){
          Some(x) => sum = x,
          None => return NTResult::Overflow,
        }
      }
    NTResult::Eval(sum)
    }

    fn mangoldt(&self) -> f64 {
      if self.reducible(){
        return (*self as u32).mangoldt()
      }
      let base = self.max_exp().0;
       if base.is_prime(){
         return (base as f64).ln()
       }
        return 0f64
    }
    
    fn mobius(&self) -> i8 {
      if self.reducible(){
        return (*self as u32).mobius()
      }
      let fctr = self.factor();
      for i in 0..fctr.len()/2{
        if fctr[2*i+1] == 2{
         return 0
        }
      }
      if fctr.len()&1 == 1{
        return -1
      }
      return 1
    }

    fn jacobi(&self, k: &Self) -> i8 {
        let mut n = *self;
        let mut p = *k;
        let mut t = 1i8;
        n %= p;

        while n != 0 {
            let zeros = n.trailing_zeros();
            n >>= zeros;

            if (p % 8 == 3 || p % 8 == 5) && (zeros % 2 == 1) {
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
        } else {
            0
        }
    }

    fn checked_jacobi(&self, k: &Self) -> NTResult<i8> {
        if k > &0 && *k % 2 == 1 {
            return NTResult::Eval(self.jacobi(k));
        }
        NTResult::Undefined
    }
    
     fn kronecker(&self, k: &Self) -> i8{
     let x = self.clone();
     if *k == 0{
      if x == 1{
         return 1
      }
     return 0
    }
   if *k == 1{
      return 1
   }
   let fctr = k.factor();
   let mut start = 0;
   let mut res = 1;
   
   if fctr[0] ==  2{
     start = 1;
     if x&1 == 0{
     res = 0;
     }
     else if x % 8 == 1 || x % 8 == 7{
      res=1
     }
     else{
       res = (-1i8).pow(fctr[1] as u32)
     }
   }
   if fctr[0] == 2 && fctr.len() == 2{
     return res
   }
   for i in start..fctr.len()/2{
     res*=self.legendre(&fctr[2*i]).pow(fctr[2*i+1] as u32);
   }
   return res
}

  fn smooth(&self) -> NTResult<Self> {
       if *self == 0{
         return NTResult::Infinite
       }
       if *self == 1{
        return NTResult::DNE
       }
        let k = self.factor();
        NTResult::Eval(k[k.len() - 2])
    }
    
    

    fn is_smooth(&self, b: &Self) -> bool {
     match self.smooth(){
      NTResult::Infinite => false,
      NTResult::Eval(x) => x <= *b, 
      _=> false,
     }
   }
}

fn delta_u64(x: u64, y: u64) -> u64 {
    //replaceable by abs_diff
    if x > y {
        x - y
    } else {
        y - x
    }
}
// mod sqr plus 1
fn mod_sqr1_64(x: u64, n: u64) -> u64 {
    ((x as u128 * x as u128 - 1) % n as u128) as u64
}
/*
// 64-bit pollard rho
fn rho_64(n: u64) -> u64 {
  let mut fact = 1u64;

  loop{
    let mut x = 2;//1152014345543559557,1103577440654186261,
    let mut y = 2;
    let mut d = 1;

    while d == 1 {
        x = mod_sqr1_64(x, n);
        y = mod_sqr1_64(mod_sqr1_64(y, n), n) % n;
        d = delta_u64(x, y).gcd(&n)
    }
    d
}
*/
fn rho_64(n: u64) -> u64 {
    // Guarantees that the result is a prime factor

    let mut x = 2; 
    let mut y = 2;
    let mut d = 1;
    loop {
        while d == 1 {
            x = mod_sqr1_64(x, n);
            y = mod_sqr1_64(mod_sqr1_64(y, n), n) % n;
            d = delta_u64(x, y).gcd(&n)
        }

        if d.is_prime() {
            return d;
        }
        d = 1; //reset for next loop
        x = u64::rng();
        y = x;
    }
}
