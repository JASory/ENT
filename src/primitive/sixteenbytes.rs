use crate::arithmetic::inlineops::*;
use crate::arithmetic::mpz::Mpz;
use crate::data::primes::DET_MAX;
use crate::data::primes::PRIMELIST;
use crate::data::primes::PRIME_INV_128;
use crate::montgomery::*;
use crate::traits::NumberTheory;
use crate::traits::Reduction;
use crate::result::NTResult;

impl NumberTheory for u128 {
    fn rng() -> Self {
        fuse(rng_64(), rng_64())
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
        sprp_128(*self, *base)
    }

    fn is_prime(&self) -> bool {
        if *self < u64::MAX as u128 {
            // truncate to faster test
            return (*self as u64).is_prime();
        }

        if self & 1 == 0 {
            return false;
        }

        if self < &DET_MAX {
            for i in PRIME_INV_128[..128].iter(){
                if (*self).wrapping_mul(*i) < *self {
                    return false;
                }
            }
            if !sprp_128(*self, 2) {
                return false;
            }
            if !sprp_128(*self, 552_491_497){
              return false;
            }
            return true;
        } else if self < &0x5A2553748E42E7B3F08195A7F78C80{
            for i in PRIME_INV_128[..128].iter() {
                if (*self).wrapping_mul(*i) < *self {
                    return false;
                }
            }

            if !sprp_128(*self, 2) {
                return false;
            }

            if detect_pseudo_128(*self) {
                return false;
            }

            if !jacobi_check_128(*self) {
                return false;
            }

            for _ in 0..10{
                if !sprp_128(*self, u128::rng() % (*self - 3) + 3) {
                    return false;
                }
            }

            return true;
        } else {
            // values greater than 2^128/727

            for i in PRIMELIST[1..64].iter() {
                if *self % *i as u128 == 0 {
                    return false;
                }
            }

            if !sprp_128(*self, 2) {
                return false;
            }

            if detect_pseudo_128(*self) {
                return false;
            }

            if !jacobi_check_128(*self) {
                return false;
            }

            for _ in 0..10 {
                if !sprp_128(*self, u128::rng() % (*self - 3) + 3) {
                    return false;
                }
            }
            /*
            if !mont_lucas(*self){
              return false
            }
              */
            return true;
        }
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

        if hi < u128::MAX {
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
        let mut count = 0u128;
        let mut start = 0u128;
        
        
        if *self < 425656284035217743{
          return (*self as u64).nth_prime().map(|y| y as u128)
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
         return (*self as u64).pi() as u128
       }
        let mut count = 0u128;
        for i in 0u128..*self {
            if i.is_prime() {
                count += 1;
            }
        }
        count
    }

    fn prime_gen(k: u32) -> NTResult<Self> {
        if k > 128 {
            return NTResult::Overflow;
        }
        if k < 65 {
            return u64::prime_gen(k).map(|x| x as u128);
        }
        let form = (1 << (k - 1)) + 1;
        let bitlength = form - 2;

        loop {
            let p = u128::rng();
            if ((p & bitlength) | form).is_prime() {
                return NTResult::Eval((p & bitlength) | form);
            }
        }
    }

    fn factor(&self) -> Vec<Self> {
        if self < &(u64::MAX as u128){
          return (*self as u64).factor().iter().map(|x| *x as u128).collect::<Vec<u128>>()
        }
        let mut n = *self;
        let twofactors = n.trailing_zeros();
        n >>= twofactors;

        let mut factors: Vec<u128> = vec![];

        if twofactors > 0 {
            factors.push(2u128);
            factors.push(twofactors as u128);
        }

        for i in PRIMELIST[1..].iter() {
            // strips out small primes
            if n % *i as u128 == 0 {
                factors.push(*i as u128);
                let mut count = 0u128;
                while n % *i as u128 == 0 {
                    count += 1;
                    n /= *i as u128;
                }
                factors.push(count);
            }
        }

        if n < u64::MAX as u128 {
            //

            let large_fact = (n as u64).factor();
            for i in large_fact {
                factors.push(i as u128)
            }
            factors
        } else {
            let large_fact = Mpz::from_u128(n).factor();

            for i in large_fact {
                factors.push(i.to_u128().unwrap())
            }
            factors
        }
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
        if *self < 0x10000000000000000 {
            return ((*self as u64).sqrt().0 as u128, 0);
        }

        let mut est = (*self as f64).sqrt() as Self + 1;

        loop {
            let s = est;
            let t = s + *self / s;
            est = t >> 1;
            println!("s: {} est: {} t: {}", s,est,t);
            if est >= s {
            if *self - (est*est) == 0{
              return (est,0)
            }
                return (s, 0);
            }
        }
    }
  // FIXME Attempt to divide by zero for x == 2^127
    fn nth_root(&self, n: &Self) -> (Self, Self) {
        if *n > 127 {
            return (1, 0);
        }

        if *n == 1 {
            return (*n, 0);
        }
        if *n == 0 {
            panic!("No integer is a zeroth factor ")
        }
        
        let mut est = (*self as f64).powf((*n as f64).recip()) as Self;
        loop {
            let s = est;
            let t = (*n - 1) * s + *self / s.pow(*n as u32 - 1);
            est = t / *n;
            if est >= s{
                return (s, 0);
            }
        }
    }

  fn max_exp(&self) -> (Self,Self){
    
      for i in 1..128{
      let p = 128-i;
      let base = self.nth_root(&p).0;
         if base.pow(p as u32) == *self{
           return(base,p)
         }
      }
     return (*self,1)
    }    

    fn radical(&self) -> NTResult<Self> {
       if self.reducible(){
         return (*self as u64).radical().map(|kishum| kishum as u128)
      }
        self.checked_factor().map(|y| y.iter().step_by(2).product::<Self>())
    }
    
    fn k_free(&self, k: &Self) -> bool {
       if self.reducible(){
         return (*self as u64).k_free(&(*k as u64))
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
        let mut gcd: u128 = *self;
        let mut new_r: u128 = *other;
        let mut bezout_1: u128 = 1;
        let mut new_s: u128 = 0;
        let mut bezout_2: u128 = 0;
        let mut new_t: u128 = 1;

        while new_r != 0 {
            let quotient = gcd / new_r;
            let mut temp: u128 = new_r;
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
      if *self < u64::MAX as u128{
        return (*self as u64).euler_totient() as u128
      }
        let factors = self.factor();
        let numerator = factors
            .iter()
            .step_by(2)
            .map(|x| x - 1u128)
            .product::<u128>();
        let denominator = factors.iter().step_by(2).product::<u128>();
        (self / denominator) * numerator
    }
    
     fn carmichael_totient(&self) -> NTResult<Self>{
     
       if  *self < u64::MAX as u128{
         return (*self as u64).carmichael_totient().map(|x| x as u128)
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

    fn jordan_totient(&self, k: &Self) -> NTResult<Self> {
        if *k > u32::MAX as u128 {
            return NTResult::Overflow;
        }
        
        if *self < 2{
           return NTResult::Eval(*self)
        }
        
        let (coef, flag) = self.overflowing_pow(*k as u32);

        if flag {
            return NTResult::CompOverflow;
        }

        let mut denom = 1u128;
        let mut numer = 1u128;

        for i in self.factor().iter().step_by(2) {
            let pow = i.pow(*k as u32);

            denom = denom * pow;

            numer *= pow - 1;
        }

        NTResult::Eval(numer * (coef / denom))
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
        if n.is_power_of_two() {
            return self.wrapping_mul(*self) & (n - 1);
        }
        pow_128(*self, 2, *n)
    }
    
   fn checked_quadratic_residue(&self, n: &Self) -> NTResult<Self> {
        if n == &0 {
            return NTResult::from_option(self.checked_mul(*self),NTResult::Overflow)
        }
        if n.is_power_of_two() {
            return NTResult::Eval(self.wrapping_mul(*self) & (n - 1));
        }
        NTResult::Eval(pow_128(*self, 2, *n))
    }

    fn product_residue(&self, other: &Self, n: &Self) -> Self {
        if n == &0 {
            return self.wrapping_mul(*other)
        }
        if n.is_power_of_two() {
            return self.wrapping_mul(*other) & (n - 1);
        }
        Mpz::from_u128(*self)
            .ref_product(&Mpz::from_u128(*other))
            .ref_euclidean(&Mpz::from_u128(*n))
            .1
            .to_u128()
            .unwrap()
    }
    
    fn checked_product_residue(&self, other: &Self, n: &Self) -> NTResult<Self>{
        if n == &0 {
            return NTResult::from_option(self.checked_mul(*self),NTResult::Overflow)
        }
        if n.is_power_of_two() {
            return NTResult::Eval(self.wrapping_mul(*other) & (n - 1));
        }
        NTResult::Eval(Mpz::from_u128(*self)
            .ref_product(&Mpz::from_u128(*other))
            .ref_euclidean(&Mpz::from_u128(*n))
            .1
            .to_u128()
            .unwrap())
    }

    fn exp_residue(&self, p: &Self, modulus: &Self) -> Self {
        if modulus == &0 {
            return self.pow(*p as u32);
        }
        pow_128(*self, *p, *modulus)
    }

    fn checked_exp_residue(&self, p: &Self, modulus: &Self) -> NTResult<Self> {
        if modulus == &0 {
            if *p > u32::MAX as u128 {
                return NTResult::Overflow;
            }
            match self.checked_pow(*p as u32) {
                Some(x) => return NTResult::Eval(x),
                None => return NTResult::Overflow,
            };
        }
        NTResult::Eval(pow_128(*self, *p, *modulus))
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
        return (*self as u64).liouville()
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
       let mut sum : u128 = 0;
       
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
        return (*self as u64).mangoldt()
      }
      let base = self.max_exp().0;
       if base.is_prime(){
         return (base as f64).ln()
       }
        return 0f64
    }
    
    fn mobius(&self) -> i8 {
    if self.reducible(){
         return (*self as u64).mobius()
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

#[inline(always)]
pub(crate) const fn split_to_u128(x: u128) -> (u128, u128) {
    (split(x).1 as u128, split(x).0 as u128)
}

#[inline] // hi,lo
pub const fn overflowing_add(lhs: (u128, u128), rhs: (u128, u128)) -> ((u128, u128), bool) {
    let (lo, carry) = lhs.1.overflowing_add(rhs.1);
    let (hi, of1) = lhs.0.overflowing_add(rhs.0);
    let (hi, of2) = hi.overflowing_add(carry as u128);
    ((hi, lo), of1 || of2)
}

pub fn u256mod128(lhs: (u128, u128), rhs: u128) -> u128 {
    if lhs.0 < rhs {
        // The result fits in 128 bits.
        div_rem1(lhs, rhs).1
    } else {
        div_rem1((lhs.0 % rhs, lhs.1), rhs).1
    }
}

pub fn u256prod(lhs: u128, rhs: u128) -> (u128, u128) {
    // hi,low
    let ((x1, x0), (y1, y0)) = (split_to_u128(lhs), split_to_u128(rhs));

    let z2 = x1 * y1;
    let (c0, z0) = split_to_u128(x0 * y0);
    let (c1, z1) = split_to_u128(x1 * y0 + c0);
    let z2 = z2 + c1;
    let (c1, z1) = split_to_u128(x0 * y1 + z1);

    (z2 + c1, z0 | z1 << 64) // hi,lo returned
}

pub fn _u256sqr(x: u128) -> (u128, u128) {
    // hi,low
    let (x1, x0) = split_to_u128(x);

    let z2 = x1 * x1;
    let m = x1 * x0;
    let (c0, z0) = split_to_u128(x0 * x0);
    let (c1, z1) = split_to_u128(m + c0);
    let z2 = z2 + c1;
    let (c1, z1) = split_to_u128(m + z1);
    (z2 + c1, z0 | z1 << 64)
}

fn _leading_zeros(pair: (u128, u128)) -> u32 {
    //hi.lo
    if pair.0 == 0 {
        pair.1.leading_zeros() + 64
    } else {
        pair.0.leading_zeros()
    }
}

fn mut_shl(pair: &mut (u128, u128), shift: u32) {
    match shift {
        0 => {}
        s if s >= 128 => {
            pair.0 = pair.1 << (s - 128);
            pair.1 = 0;
        }
        s => {
            pair.0 <<= s;
            pair.0 |= pair.1 >> (128 - s);
            pair.1 <<= s;
        }
    }
}
// which is all that is necessary until a faster optimization is constructed.
// Any library that uses this function as a general euclidean function is to be considered critically broken.
pub fn div_rem1(pair: (u128, u128), other: u128) -> (u128, u128) {
    //takes hi,lo pair

    const RADIX: u128 = 0x10000000000000000;

    // Normalize the divisor.
    let s = other.leading_zeros();
    let d = other << s; // numerator, denominator
    let mut zqk = pair;
    mut_shl(&mut zqk, s);
    let p = zqk;
    let (d1, d0) = split_to_u128(d);
    let (n1, n0) = split_to_u128(p.1); // split lower part of dividend

    let (mut q1, mut rhat) = p.0.euclidean_div(&d1);

    while q1 >= RADIX || q1 * d0 > RADIX * rhat + n1 {
        q1 -= 1;
        rhat += d1;
        if rhat >= RADIX {
            break;
        }
    }

    let r21 =
        p.0.wrapping_mul(RADIX)
            .wrapping_add(n1)
            .wrapping_sub(q1.wrapping_mul(d));

    // Compute the second quotient digit q0.
    let (mut q0, mut rhat) = r21.euclidean_div(&d1);

    // q0 has at most error 2. No more than 2 iterations.
    while q0 >= RADIX || q0 * d0 > RADIX * rhat + n0 {
        q0 -= 1;
        rhat += d1;
        if rhat >= RADIX {
            break;
        }
    }

    let r = (r21
        .wrapping_mul(RADIX)
        .wrapping_add(n0)
        .wrapping_sub(q0.wrapping_mul(d)))
        >> s;
    let q = q1 * RADIX + q0;
    (q, r)
}

/*  Replacement function
pub fn rem(pair: (u128,u128), other: u128) -> u128 {

   const RADIX : u128 = 0x10000000000000000;

   // Normalize the divisor.
    let s = other.leading_zeros();
    let n = other << s; // numerator, denominator
    let mut zqk = pair;
    mut_shl(&mut zqk, s);

    let (n_hi, nlo) = split_to_u128(n);

    let (hi_hi, hi_lo) = split(pair.0);

    let (lo_hi, lo_lo) = split(pair.1);

    let quo =  divide3by2(hi_hi, hi_lo, lo_hi, n_hi, n_lo);


}

*/

fn detect_pseudo_128(x: u128) -> bool {
    //
    for i in 2..32 {
        let fctr = 2 * i;
        let sq = (x - 1) / fctr;
        let k = sq.sqrt().0;
        if ((k * k + k) * fctr + k + 1) == x {
            return true;
        }
    }
    return false;
}

fn jacobi_check_128(x: u128) -> bool {
    let mut witness = 3u128;
    loop {
        if witness.jacobi(&x) == -1 {
            break;
        }
        witness += 1;
    }
    x.is_sprp(&witness)
}


