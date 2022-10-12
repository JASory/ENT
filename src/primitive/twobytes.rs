use crate::data::nt_data::LIOUVILLE_LUT;
use crate::data::primes::PRIMELIST;
use crate::traits::NumberTheory;

use crate::arithmetic::inlineops::*;

impl NumberTheory for u16 {
    fn rng() -> Self {
        (rng_32() >> 16) as u16
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
        if *self < u8::MAX as u16 {
            // tree down to u8 if it fits
            return (*self as u8).is_prime();
        }

        for i in PRIMELIST[..54].iter() {
            if self == i {
                return true;
            }
            if self % i == 0 {
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
                if witness.euclid_gcd(self) == 1 {
                    break 'witness;
                }
                witness += 1;
            }

            if witness.exp_residue(&x_minus, self) != 1 {
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

        if hi < u16::MAX {
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

    fn nth_prime(&self) -> Option<Self> {
        let mut count = 0u16;
        let mut start = 0u16;
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
        let mut count = 0u16;
        for i in 0u16..*self {
            if i.is_prime() {
                count += 1;
            }
        }
        count
    }

    fn prime_gen(k: u32) -> u16 {
        if k > 16 {
            panic!("Outside the limit of the datatype")
        }
        if k < 8 {
            return u8::prime_gen(k) as u16;
        }
        let form = (1 << (k - 1)) + 1;
        let bitlength = form - 2;

        loop {
            let q = u64::rng();
            let p = unsafe { std::mem::transmute::<u64, (u16, u16, u16, u16)>(q) };

            if ((p.0 & bitlength) | form).is_prime() {
                return (p.0 & bitlength) | form;
            }
            if ((p.1 & bitlength) | form).is_prime() {
                return (p.1 & bitlength) | form;
            }
            if ((p.2 & bitlength) | form).is_prime() {
                return (p.2 & bitlength) | form;
            }
            if ((p.3 & bitlength) | form).is_prime() {
                return (p.3 & bitlength) | form;
            }
        }
    }

    fn factor(&self) -> Vec<Self> {
    
      if self < &255{
        return (*self as u8).factor().iter().map(|x| *x as u16).collect::<Vec<u16>>()
      }
        let mut n = *self;
        let twofactors = n.trailing_zeros();
        n >>= twofactors;

        let mut factors: Vec<u16> = vec![];

        if twofactors > 0 {
            factors.push(2);
            factors.push(twofactors as u16);
        }

        for i in PRIMELIST[1..54].iter() {
            // strips out small primes
            if n % i == 0 {
                factors.push(*i);
                let mut count = 0u16;
                while n % i == 0 {
                    count += 1;
                    n /= i;
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

    fn sqrt(&self) -> (Self, Self) {
        ((*self as f64).sqrt() as Self, 0)
    }

    fn nth_root(&self, n: &Self) -> (Self, Self) {
        if *n > 15 {
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

    fn radical(&self) -> Self {
        self.factor().iter().step_by(2).product::<u16>()
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
        let mut gcd: u16 = *self;
        let mut new_r: u16 = *other;
        let mut bezout_1: u16 = 1;
        let mut new_s: u16 = 0;
        let mut bezout_2: u16 = 0;
        let mut new_t: u16 = 1;

        while new_r != 0 {
            let quotient = gcd / new_r;
            let mut temp: u16 = new_r;
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
        (*self / cf) * (*other)
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
        let numerator = factors.iter().step_by(2).map(|x| x - 1u16).product::<u16>();
        let denominator = factors.iter().step_by(2).product::<u16>();
        (self / denominator) * numerator
    }

    fn jordan_totient(&self, k: &Self) -> Option<Self> {
        let (coef, flag) = self.overflowing_pow(*k as u32);

        if flag {
            return None;
        }

        let mut denom = 1u16;
        let mut numer = 1u16;

        for i in self.factor().iter().step_by(2) {
            let pow = i.pow(*k as u32);

            denom *= pow;

            numer *= pow - 1;
        }

        Some(numer * (coef / denom))
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
            match self.checked_mul(*self) {
                Some(x) => return x,
                None => panic!("Element of residue class exceeds datatype bound"),
            };
        }
        ((*self as u32 * *self as u32) % *n as u32) as u16
    }

    fn product_residue(&self, other: &Self, n: &Self) -> Self {
        if n == &0 {
            match self.checked_mul(*other) {
                Some(x) => return x,
                None => panic!("Element of residue class exceeds datatype bound"),
            };
        }
        ((*self as u32 * *other as u32) % *n as u32) as u16
    }

    fn exp_residue(&self, p: &Self, modulus: &Self) -> Self {
        if modulus == &0 {
            match self.checked_pow(*p as u32) {
                Some(x) => return x,
                None => panic!("Element of residue class exceeds datatype bound"),
            };
        }

        if modulus.is_power_of_two() {
            return self.wrapping_pow((*p) as u32) & (*modulus - 1);
        }

        let mut z = 1u32;
        let mut base = *self as u32;
        let n = *modulus as u32;
        let mut pow = *p;
        if pow == 0 {
            return z as u16;
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

        (base * z % n) as u16
    }

    fn checked_exp_residue(&self, p: &Self, modulus: &Self) -> Option<Self> {
        if modulus == &0 {
            match self.checked_pow(*p as u32) {
                Some(x) => return Some(x),
                None => return None,
            };
        }

        if modulus.is_power_of_two() {
            return Some(self.wrapping_pow((*p) as u32) & (*modulus - 1));
        }

        let mut z = 1u32;
        let mut base = *self as u32;
        let n = *modulus as u32;
        let mut pow = *p;
        if pow == 0 {
            return Some(z as u16);
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

        Some((base * z % n) as u16)
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
        if (LIOUVILLE_LUT[(*self / 64) as usize] >> (*self % 64)) & 1 == 1 {
            return -1;
        }
        return 1;
    }
    
    fn derivative(&self) -> Option<Self> {
       let fctr = self.factor();
       let mut sum : u16 = 0;
       
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
