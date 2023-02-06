use crate::data::primes::PRIMELIST;
use crate::data::primes::PRIME_INV_64;
use crate::traits::NumberTheory;

use crate::arithmetic::inlineops::*;
use crate::data::hashtable::BASE_32;
use crate::montgomery::*;
use crate::sieve::*;

impl NumberTheory for u32 {
    fn rng() -> Self {
        rng_32()
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
        sprp_32(*self, *base)
    }

    fn is_prime(&self) -> bool {
        if *self < u16::MAX as u32 {
            // tree down to u16 if it fits
            return (*self as u16).is_prime();
        }

        if *self & 1 == 0 {
            return false;
        }

        for i in PRIME_INV_64[..128].iter() {
            if (*self as u64).wrapping_mul(*i) < *self as u64 {
                return false;
            }
        }

        let idx = (self.wrapping_mul(0xA9DDE81F) / 0x147AE15) as usize;

        self.is_sprp(&(BASE_32[idx] as u32))
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
                if witness.euclid_gcd(&self) == 1 {
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
        let hi = std::cmp::max(*self, *sup);

        prime_list_32(inf as usize, hi as usize, &erasto_sieve(hi as usize))
    }

    fn nth_prime(&self) -> Option<Self> {
        let mut count = 0u32;
        let mut start = 0u32;
        loop {
            start += 1;

            if start == Self::MAX {
                return None;
            }
            if start.is_prime() {
                count += 1;
            }
            if count == *self {
                return Some(start);
            }
        }
    }

    fn pi(&self) -> Self {
       if *self < 9{
         return (*self as u8).pi() as u32
       }
        prime_count((*self) as usize, &erasto_sieve(*self as usize))
    }

    fn prime_gen(k: u32) -> Option<Self> {
        if k > 32 {
            return None;
        }
        if k < 16 {
            return u16::prime_gen(k).map(|x| x as u32);
        }
        let form = (1 << (k - 1)) + 1;
        let bitlength = form - 2;

        loop {
            let q = u64::rng();
            let p = unsafe { std::mem::transmute::<u64, (u32, u32)>(q) };

            if ((p.0 & bitlength) | form).is_prime() {
                return Some((p.0 & bitlength) | form);
            }
            if ((p.1 & bitlength) | form).is_prime() {
                return Some((p.1 & bitlength) | form);
            }
        }
    }

    fn factor(&self) -> Vec<Self> {
      if self < &65535{
      return (*self as u16).factor().iter().map(|x| *x as u32).collect::<Vec<u32>>()
      }
        let mut n = *self;
        let twofactors = n.trailing_zeros();
        n >>= twofactors;

        let mut factors: Vec<u32> = vec![];

        if twofactors > 0 {
            factors.push(2);
            factors.push(twofactors);
        }

        if n.is_prime() {
            factors.push(n);
            factors.push(1);
            return factors;
        }

        for i in PRIMELIST[1..].iter() {
            // strips out small primes
            if n % *i as u32 == 0 {
                factors.push(*i as u32);
                let mut count = 0u32;
                while n % *i as u32 == 0 {
                    count += 1;
                    n /= *i as u32;
                }
                factors.push(count);
            }
        }

        if n == 1 {
            return factors;
        }

        while n != 1 {
            let k = rho_32(n);
            factors.push(k);
            let mut count = 0u32;
            while n % k == 0 {
                count += 1;
                n /= k;
            }
            factors.push(count);
        }
        factors
    }
    
    fn checked_factor(&self) -> Option<Vec<Self>>{
     
     if *self < 2{
       return None
     }
     
     Some(self.factor())
    
    }


    fn sqrt(&self) -> (Self, Self) {
        ((*self as f64).sqrt() as Self, 0)
    }

    fn nth_root(&self, n: &Self) -> (Self, Self) {
        if *n > 31 {
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

    
    fn radical(&self) -> Option<Self> {
        self.checked_factor().map(|y| y.iter().step_by(2).product::<Self>())
    }

    fn k_free(&self, k: &Self) -> bool {
        let factors = self.factor();
        for (idx, el) in factors.iter().enumerate() {
            if el == k && idx != 0 {
                return false;
            }
        }
        true
    }

    fn euclid_gcd(&self, other: &Self) -> Self {
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
        let mut gcd: u32 = *self;
        let mut new_r: u32 = *other;
        let mut bezout_1: u32 = 1;
        let mut new_s: u32 = 0;
        let mut bezout_2: u32 = 0;
        let mut new_t: u32 = 1;

        while new_r != 0 {
            let quotient = gcd / new_r;
            let mut temp: u32 = new_r;
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
        let cf = self.euclid_gcd(other);
        (*self / cf) * *other
    }

    fn checked_lcm(&self, other: &Self) -> Option<Self> {
         if self == &0 && other == &0{
           return Some(0)
         }
        let cf = self.euclid_gcd(other);
        let (v, flag) = (*self / cf).overflowing_mul(*other);
        if flag {
            return None;
        }
        Some(v)
    }

    fn euler_totient(&self) -> Self {
        let factors = self.factor();
        let numerator = factors.iter().step_by(2).map(|x| x - 1u32).product::<u32>();
        let denominator = factors.iter().step_by(2).product::<u32>();
        (self / denominator) * numerator
    }

    fn jordan_totient(&self, k: &Self) -> Option<Self> {
        let (coef, flag) = self.overflowing_pow(*k);

        if flag {
            return None;
        }

        let mut denom = 1u32;
        let mut numer = 1u32;

        for i in self.factor().iter().step_by(2) {
            let pow = i.pow(*k);

            denom = denom * pow;

            numer *= pow - 1;
        }

        Some(numer * (coef / denom))
    }
    
    fn carmichael_totient(&self) -> Option<Self>{
    
       if  *self < 65535{
        return (*self as u16).carmichael_totient().map(|x| x as u32)
       }
       
       let fctr = self.factor();
       let base = fctr.iter().step_by(2).map(|z| *z).collect::<Vec<Self>>();
       let mut result = 1;
      for (idx,el) in base.iter().enumerate(){
        if el == &2 && fctr[1] > 2{
         let phi = ((el.pow(fctr[2*idx+1]) /el) *(el-1)) /2;
          result = result.lcm(&phi);
        }
       else{
         let phi =  (el.pow(fctr[2*idx+1])/el)*(el-1);
         result = result.lcm(&phi);
       } 
      }
     Some(result)
    }

    fn dedekind_psi(&self, k: &Self) -> Option<Self> {
        let (k2, flag) = k.overflowing_shl(1);
        if flag {
            return None;
        }
        match self.jordan_totient(&k2) {
            Some(y) => Some(y / self.jordan_totient(k).unwrap()),
            None => None,
        }
    }

   fn quadratic_residue(&self, n: &Self) -> Self {
        if n == &0 {
            return self.wrapping_mul(*self)
        }
        ((*self as u64 * *self as u64) % *n as u64) as Self
    }
    
    fn checked_quadratic_residue(&self, n: &Self) -> Option<Self> {
        if n == &0 {
            return self.checked_mul(*self)
        }
        Some(((*self as u64 * *self as u64) % *n as u64) as Self)
    }
    

    fn product_residue(&self, other: &Self, n: &Self) -> Self {
        if n == &0 {
            return self.wrapping_mul(*other)
        }
        ((*self as u64 * *other as u64) % *n as u64) as Self
    }
    
    fn checked_product_residue(&self, other: &Self, n: &Self) -> Option<Self> {
        if n == &0 {
            return self.checked_mul(*other)
        }
        Some(((*self as u64 * *other as u64) % *n as u64) as Self)
    }

    fn exp_residue(&self, p: &Self, modulus: &Self) -> Self {
        if modulus == &0 {
          return self.wrapping_pow(*p)
        }

        if modulus.is_power_of_two() {
            return self.wrapping_pow(*p) & (*modulus - 1);
        }

        pow_32(*self, *p, *modulus)
    }

    fn checked_exp_residue(&self, p: &Self, modulus: &Self) -> Option<Self> {
        if modulus == &0 {
            match self.checked_pow(*p) {
                Some(x) => return Some(x),
                None => return None,
            };
        }

        if modulus.is_power_of_two() {
            return Some(self.wrapping_pow(*p) & (*modulus - 1));
        }

        Some(pow_32(*self, *p, *modulus))
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

    fn checked_legendre(&self, p: &Self) -> Option<i8> {
        if p == &2 || !p.is_prime() {
            return None;
        }
        Some(self.legendre(p))
    }

    fn liouville(&self) -> i8 {
        let primeomega = self.factor()[1..].iter().step_by(2).sum::<Self>();
        if primeomega & 1 == 0 {
            return 1;
        }
        return -1;
    }
    
    fn derivative(&self) -> Option<Self> {
       let fctr = self.factor();
       let mut sum : u32 = 0;
       
     for i in 0..fctr.len() / 2 {
        match sum.checked_add(fctr[2 * i + 1] * (*self / fctr[2 * i])){
          Some(x) => sum = x,
          None => return None,
        }
      }
    Some(sum)
    }

    fn mangoldt(&self) -> f64 {
        let fctr = self.factor();
        if fctr.len() != 2 {
            return 0f64;
        }
        (fctr[0] as f64).ln()
    }
    
    fn mobius(&self) -> i8 {
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

    fn checked_jacobi(&self, k: &Self) -> Option<i8> {
        if k > &0 && *k % 2 == 1 {
            return Some(self.jacobi(k));
        }
        None
    }

    fn smooth(&self) -> Self {
        let k = self.factor();
        k[k.len() - 2]
    }

    fn is_smooth(&self, b: &Self) -> bool {
        &self.smooth() <= b
    }
}

// 32-bit pollard rho
fn delta_u32(x: u32, y: u32) -> u32 {
    // replacable by abs_diff
    if x > y {
        x - y
    } else {
        y - x
    }
}

// mod sqr plus 1
fn mod_sqr1_32(x: u32, n: u32) -> u32 {
    ((x as u64 * x as u64 + 1) % n as u64) as u32
}

// 32-bit pollard rho
fn rho_32(n: u32) -> u32 {
    let mut x = 2;
    let mut y = 2;
    let mut d = 1;

    while d == 1 {
        x = mod_sqr1_32(x, n);
        y = mod_sqr1_32(mod_sqr1_32(y, n), n) % n;
        d = delta_u32(x, y).euclid_gcd(&n)
    }
    d
}
