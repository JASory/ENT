use crate::data::nt_data::LIOUVILLE_LUT;
use crate::data::nt_data::EULER_TOTIENT_LUT;
use crate::data::nt_data::CARMICHAEL_LUT;
use crate::data::nt_data::DERIVATIVE_LUT;

use crate::data::primes::PRIMELIST;
use crate::traits::NumberTheory;

use crate::result::NTResult;

use crate::arithmetic::inlineops::*;

impl NumberTheory for u8 {
    fn rng() -> Self {
        (rng_32() >> 24) as u8
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
        let pminus = *self - 1;
        let zeroes = pminus.trailing_zeros(); // Breaks number down to p= d*2^n -1
        let d = pminus >> zeroes;
        let mut x = base.exp_residue(&d, self);

        if x == 1 || x == pminus {
            return true;
        }

        for _ in 1..zeroes {
            x = x.quadratic_residue(self);
            if x == pminus {
                return true;
            }
        }
        false // otherwise it fails
    }

    fn is_prime(&self) -> bool {
        if *self == 1 || *self == 0 {
            return false;
        }
        for i in PRIMELIST[..6].iter() {
            if *self == *i as u8 {
                return true;
            }
            if self % *i as u8 == 0 {
                return false;
            }
        }
        true
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

        if hi < u8::MAX {
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
        if *self == 0{
          return NTResult::DNE;
        }
        if *self > 54{
         return NTResult::Overflow
        }
        NTResult::Eval(PRIMELIST[*self as usize] as u8)
    }

    fn pi(&self) -> Self {
        let mut count = 0u8;
        if *self == 0{
          return 0
        }
        for i in 0u8..*self {
            if i.is_prime() {
                count += 1;
            }
        }
        count
    }

    fn prime_gen(k: u32) -> NTResult<Self> {
        if k > 8 {
            return NTResult::Overflow;
        }
        if k < 1 {
            return NTResult::DNE;
        }
        if k == 1{
          return NTResult::Eval(2)
        }

        let form = (1 << (k - 1)) + 1;
        let bitlength = form - 2;
        loop {
            // generate a list of bits then split into datatype size, rarely requires more than one iteration
            let q = u32::rng();
            let p = q.to_be_bytes();

            if ((p[0] & bitlength) | form).is_prime() {
                return NTResult::Eval((p[0] & bitlength) | form);
            }
            if ((p[1] & bitlength) | form).is_prime() {
                return NTResult::Eval((p[1] & bitlength) | form);
            }
            if ((p[2] & bitlength) | form).is_prime() {
                return NTResult::Eval((p[2] & bitlength) | form);
            }
            if ((p[3] & bitlength) | form).is_prime() {
                return NTResult::Eval((p[3] & bitlength) | form);
            }
        }
    }

    fn factor(&self) -> Vec<u8> {

        let mut n = *self;
        let twofactors = n.trailing_zeros();
        n >>= twofactors;

        let mut factors: Vec<u8> = vec![];

        if twofactors > 0 {
            factors.push(2);
            factors.push(twofactors as u8);
        }

        for i in PRIMELIST[1..6].iter() {
            // strips out small primes
            if n % *i as u8 == 0 {
                factors.push(*i as u8);
                let mut count = 0u8;
                while n % *i as u8 == 0 {
                    count += 1;
                    n /= *i as u8;
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
        ((*self as f64).sqrt() as Self, 0)
    }

    fn nth_root(&self, n: &Self) -> (Self, Self) {
        if *n > 7 {
            return (1, 0);
        }

        if *n == 1 {
            return (*n, 0);
        }

        if *n == 0 {
            panic!("No integer is a zeroth factor ")
        }
        (((*self as f64).powf((*n as f64).recip())) as Self, 0)
    }
    
    fn max_exp(&self) -> (Self,Self){
    
      for i in 1..8{
      let p = 8-i;
      let base = self.nth_root(&p).0;
         if base.pow(p as u32) == *self{
           return(base,p)
         }
      }
     return (*self,1)
    }    

    
    fn radical(&self) -> NTResult<Self> {
       if *self == 0{
          return NTResult::Infinite
       }
       if *self == 1{
          return NTResult::Eval(1)
       }
        self.checked_factor().map(|y| y.iter().step_by(2).product::<Self>())
    }

    fn k_free(&self, k: &Self) -> bool {
        if *self == 0{
          return false 
        }
        if *self == 1{
         return true
        }
        let factors = self.factor();
        for (idx, fact) in factors.iter().enumerate() {
            if fact == k && idx > 0 {
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
        let mut gcd: u8 = *self;
        let mut new_r: u8 = *other;
        let mut bezout_1: u8 = 1;
        let mut new_s: u8 = 0;
        let mut bezout_2: u8 = 0;
        let mut new_t: u8 = 1;

        while new_r != 0 {
            let quotient = gcd / new_r;
            let mut temp: u8 = new_r;
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

    fn euler_totient(&self) -> u8 {
        EULER_TOTIENT_LUT[*self as usize]
    }

    fn jordan_totient(&self, k: &Self) -> NTResult<Self> {
        if *self < 2{
           return NTResult::Eval(*self)
        }
        let (coef, flag) = self.overflowing_pow(*k as u32);
        
        if flag {
            return NTResult::CompOverflow;
        }

        let mut denom = 1u8;
        let mut numer = 1u8;

        for i in self.factor().iter().step_by(2) {
            let pow = i.pow(*k as u32);

            denom = denom * pow;

            numer *= pow - 1;
        }

        NTResult::Eval(numer * (coef / denom))
    }
    
    fn carmichael_totient(&self) -> NTResult<Self>{
        if *self == 0{
         return NTResult::Infinite;
        }
        NTResult::Eval(CARMICHAEL_LUT[*self as usize])
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
        ((*self as u16 * *self as u16) % *n as u16) as Self
    }
    
    fn checked_quadratic_residue(&self, n: &Self) -> NTResult<Self> {
        if n == &0 {
            return NTResult::from_option(self.checked_mul(*self),NTResult::Overflow)
        }
        NTResult::Eval( ((*self as u16 * *self as u16) % *n as u16) as Self)
    }

    fn product_residue(&self, other: &Self, n: &Self) -> Self {
        if n == &0 {
            return self.wrapping_mul(*other)
        }
        ((*self as u16 * *other as u16) % *n as u16) as Self
    }
    
    fn checked_product_residue(&self, other: &Self, n: &Self) -> NTResult<Self> {
        if n == &0 {
            return NTResult::from_option(self.checked_mul(*self),NTResult::Overflow)
        }
        NTResult::Eval( ((*self as u16 * *other as u16) % *n as u16) as Self) 
    }

    fn exp_residue(&self, p: &Self, modulus: &Self) -> Self {
        if modulus == &0 {
           return (*self).wrapping_pow(*p as u32)
        }
        if modulus.is_power_of_two() {
            return self.wrapping_pow((*p) as u32) & (*modulus - 1);
        }
        let mut z = 1u16;
        let mut base = *self as u16;
        let n = *modulus as u16;
        let mut pow = *p;
        if pow == 0 {
            return z as u8;
        }

        while pow > 1 {
            if pow % 2 == 0 {
                base = base * base % n;
                pow >>= 1;
            } else {
                z = base * z % n;
                base = base * base % n;
                pow = (pow - 1) >> 1;
            }
        }

        (base * z % n) as u8
    }

    fn checked_exp_residue(&self, p: &Self, modulus: &Self) -> NTResult<Self> {
        if modulus == &0 {
            match self.checked_pow(*p as u32) {
                Some(x) => return NTResult::Eval(x),
                None => return NTResult::Overflow,
            };
        }
        if modulus.is_power_of_two() {
            return NTResult::Eval(self.wrapping_pow((*p) as u32) & (*modulus - 1));
        }
        let mut z = 1u16;
        let mut base = *self as u16;
        let n = *modulus as u16;
        let mut pow = *p;
        if pow == 0 {
            return NTResult::Eval(z as u8);
        }

        while pow > 1 {
            if pow % 2 == 0 {
                base = base * base % n;
                pow >>= 1;
            } else {
                z = base * z % n;
                base = base * base % n;
                pow = (pow - 1) >> 1;
            }
        }

        NTResult::Eval((base * z % n) as u8)
    }

    fn legendre(&self, p: &Self) -> i8 {
        let k = self.exp_residue(&((*p - 1) >> 1), p);
        if k == 1 {
            return 1;
        };
        if k == *p - 1 {
            return -1;
        };
        0
    }

    fn checked_legendre(&self, p: &Self) -> NTResult<i8> {
        if p == &2 || !p.is_prime() {
            return NTResult::Undefined;
        }
        NTResult::Eval(self.legendre(p))
    }
    
    fn derivative(&self) -> NTResult<Self> {
      if *self > 94{
         return NTResult::Overflow;
      }
      NTResult::Eval(DERIVATIVE_LUT[*self as usize])
    }

    fn liouville(&self) -> i8 {
        if (LIOUVILLE_LUT[(*self / 64) as usize] >> (*self % 64)) & 1 == 1 {
            return -1;
        }
        return 1;
    }
    
    fn mobius(&self) -> i8 {
      if *self == 0{
        return 0i8
      }
      if *self == 1{
        return 1
      }
      let fctr = self.factor();
      if fctr.len() == 1{ // if only one factor then return -1
         return -1
      }
      for i in 0..fctr.len()/2{
        if fctr[2*i+1]  > 1{
         return 0
        }
      }
      let fctrsum = fctr[1..].iter().step_by(2).sum::<Self>();
      if fctrsum&1 == 1{// if odd number of factors and square free
        return -1
      }
      return 1
    }

    fn mangoldt(&self) -> f64 {
      let base = self.max_exp().0;
       if base.is_prime(){
         return (base as f64).ln()
       }
        return 0f64
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

    fn checked_jacobi(&self, k: &Self) -> NTResult<i8> {
        if k > &0 && *k % 2 == 1 {
            return NTResult::Eval(self.jacobi(k));
        }
        NTResult::Undefined
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
