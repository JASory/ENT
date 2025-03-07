use crate::structs::{Certificate,Factorization};
use crate::ntrait::{NumberTheory,MontArith,Reduction};
use crate::primitive::factorprim::{poly_factor,drbg};
use crate::data::primes::PRIMELIST;
use crate::result::NTResult;

impl NumberTheory for u64 {

    fn is_unit(&self) -> bool{
	if *self == 1{
	    return true;
	}
	false
    }

    fn rng() -> Self {
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if is_x86_feature_detected!("rdrand"){// USE RDRAND if possible
            let mut x: u64 = 0;
            unsafe { core::arch::x86_64::_rdrand64_step(&mut x) };
            return x;
        }// If processor is x86 and does not support RDRAND use xor shift
       return drbg(std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos() as u64)
    }
        #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
    {// All other architectures use xor shift
        drbg(std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos() as u64)
    }
    }

    fn residue(&self, ring: Self) -> Self{
      if ring == 0{
        return *self
      }
        *self % ring
    }
     
    fn euclidean_div(&self, other: Self) -> (Self, Self) {
        (*self / other, *self % other)
    }
    
    fn mul_inverse(&self, n: Self) -> NTResult<Self>{
       let (gcd,inv,_) = self.extended_gcd(n);
       
       if gcd != 1{
          return NTResult::DNE
       }
       NTResult::Eval(inv)
    }
    
    fn fermat(&self,base: Self) -> bool{
        if *self&1 == 1{
           let one = self.n_identity();
           let mbase = base.to_mont(*self);
           let inv = self.inv_2();
           return mbase.mont_pow(one,*self-1,inv,*self)==one;
       }
       base.exp_residue(*self-1,*self)==1
       
    }

    fn strong_fermat(&self, base: Self) -> bool {
        if *self&1==0{
            return self.fermat(base);
        }
        let inv = self.inv_2();//machine_prime::mul_inv2(x);
        let tzc = self.wrapping_sub(1).trailing_zeros();
        let one = self.n_identity();
        let oneinv = machine_prime::mont_prod(machine_prime::mont_sub(*self, one, *self), one, *self, inv);
        let b = machine_prime::to_mont(base,*self);
    
        machine_prime::strong_fermat(*self, tzc, b, one, oneinv, inv)
    }

    fn is_prime(&self) -> bool {
        machine_prime::is_prime(*self)
    }

    fn prime_proof(&self) -> Certificate<Self> {

        if *self == 2 {
            return Certificate::new(*self,3,vec![]);
        }

        let x_minus = *self - 1;
        let fctrs = x_minus.factor().unwrap()
            .factor_iter().map(|x| x_minus/ *x)
            .collect::<Vec<Self>>();

        loop {
            // loops until we get a 

            let mut witness = Self::rng() % (*self - 2) + 2;
            
            'witness: loop {

                if witness.coprime(*self) {

		  break 'witness;
		}
              
                witness += 1;
            }
              
	   if !self.strong_fermat(witness){
		return Certificate::new(*self,witness,fctrs)	
	   }

           
            'inner: for (idx, i) in fctrs.iter().enumerate() {
                if witness.exp_residue(*i, *self).is_unit(){
                    break 'inner;
                }
                if idx == fctrs.len() - 1 {
                    return Certificate::new(*self,witness,fctrs);
                }
            }
	  }
    }

    fn prime_list(&self, sup: Self) -> Vec<Self> {
       
        let inf = std::cmp::min(*self, sup);
        let mut hi = std::cmp::max(*self, sup);

        hi = hi.saturating_add(1);

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
        
/*        if *self < 203280221{
          return (*self as u32).nth_prime().map(|y| y as u64)
        }
  */      
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
       
        let mut count = 203280221u64;
        for i in (1u64<<32)..*self {
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

    fn factor(&self) -> NTResult<Factorization<Self>> {
       /* if self < &4294967295{
          return (*self as u32).factor().iter().map(|x| *x as u64).collect::<Vec<u64>>()
        }
       */ 
        if *self == 0{
           return NTResult::InfiniteSet;
        }
        
        if *self == 1{
           return NTResult::DNE
        }
        let mut n = *self;
        let twofactors = n.trailing_zeros();
        n >>= twofactors;

        let mut factors: Factorization<u64> = Factorization::new();

        if twofactors > 0 {
            factors.add_factor(2);
            factors.add_power(twofactors as u64);
        }

        for i in PRIMELIST[1..128].iter() {
            // strips out small primes
            if n % *i as u64 == 0 {
                factors.add_factor(*i as u64);
                let mut count = 0u64;
                while n % *i as u64 == 0 {
                    count += 1;
                    n /= *i as u64;
                }
                factors.add_power(count);
            }
        }

        if n == 1 {
            return  NTResult::Eval(factors);
        }

        if n.is_prime() {
            factors.add_factor(n);
            factors.add_power(1);
            return  NTResult::Eval(factors);
        }

        while n != 1 {
            let k = poly_factor(n);
            factors.add_factor(k);
            let mut count = 0u64;
            while n % k == 0 {
                count += 1;
                n /= k;
            }
            factors.add_power(count);
            
            if n.is_prime() {
            factors.add_factor(n);
            factors.add_power(1);
            return  NTResult::Eval(factors);
          }
          
        }
        NTResult::Eval(factors)
    }
    


    fn sqrt(&self) -> (Self, Self) {
       // if *self < 0x100000000 {
       //     return ((*self as u32).sqrt().0 as u64, 0);
       // }
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

    fn nth_root(&self, n: Self) -> (Self, Self) {
        if n > 63 {
            return (1, 0);
        }

        if n == 1 {
            return (n, 0);
        }

        if n == 0 {
            panic!("No integer is a zeroth factor ")
        }

        let mut est = (*self as f64).powf((n as f64).recip()) as Self + 1;

        loop {
            let s = est;
            let t = (n - 1) * s + *self / s.pow(n as u32 - 1);
            est = t / n;
            if est >= s {
                return (s, 0);
            }
        }
    }

    fn max_exp(&self) -> (Self,Self){
    
      for i in 1..64{
      let p = 64-i;
      let base = self.nth_root(p).0;
         if base.pow(p as u32) == *self{
           return(base,p)
         }
      }
     (*self,1)
    }    
    
    fn radical(&self) -> NTResult<Self> {
      if self.reducible(){
         return (*self as u32).radical().map(|kishum| kishum as u64)
      }
        self.factor().map(|y| y.factor_iter().product::<Self>())
    }

    fn k_free(&self, k: Self) -> bool {
       if self.reducible(){
         return (*self as u32).k_free(k as u32)
      }
      // FIXME Remove unwrap
        let factors = self.factor().unwrap();
        factors.k_free(k)
    }

    fn gcd(&self, other: Self) -> Self {
        let mut a = *self;
        let mut b = other;
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

    fn extended_gcd(&self, other: Self) -> (Self, Self, Self) {
        let mut gcd: u64 = *self;
        let mut new_r: u64 = other;
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
            if bezout_1 < quotient.product_residue(temp, other) {
                // First bezout coefficient is computed over Z[q]
                new_s = other - (quotient.product_residue(temp, other) - bezout_1)
            } else {
                new_s = bezout_1.wrapping_sub(quotient * temp);
            }

            bezout_1 = temp;

            temp = new_t;
            if bezout_2 < quotient.product_residue(temp, *self) {
                // Second bezout coefficient is computed over Z[p]
                new_t = *self - (quotient.product_residue(temp, *self) - bezout_2)
            } else {
                new_t = bezout_2.wrapping_sub(quotient.product_residue(temp, *self));
            }
            bezout_2 = temp
        }
        (gcd, bezout_1, bezout_2)
    }
/*
    fn coprime(&self, other: Self) -> bool{
       self.gcd(other)==1
    }
*/
    fn lcm(&self, other: Self) -> NTResult<Self> {
        if self == &0 && other == 0{
           return NTResult::Eval(0)
         }
        let cf = self.gcd(other);
        let (v, flag) = (*self / cf).overflowing_mul(other);
        if flag {
            return NTResult::Overflow;
        }
        NTResult::Eval(v)
    }

    fn euler_totient(&self) -> Self {
      if self.reducible(){
        return (*self as u32).euler_totient() as u64
      }
        let factors = self.factor().unwrap();
        let numerator = factors.factor_iter().map(|x| x - 1u64).product::<u64>();
        let denominator = factors.factor_iter().product::<u64>();
        (self / denominator) * numerator
    }

    fn jordan_totient(&self, k: Self) -> NTResult<Self> {
        if k > u32::MAX as u64 {
            return NTResult::Overflow;
        }
        
        if *self < 2{
           return NTResult::Eval(*self)
        }
        
        let (coef, flag) = self.overflowing_pow(k as u32);

        if flag {
            return NTResult::CompOverflow;
        }

        let mut denom = 1u64;
        let mut numer = 1u64;

        for i in self.factor().unwrap().factor_iter(){
            let pow = i.pow(k as u32);

            denom *= pow;

            numer *= pow - 1;
        }

        NTResult::Eval(numer * (coef / denom))
    }
    
    fn exponent(&self) -> NTResult<Self>{ 
    
       if self.reducible(){
         return (*self as u32).exponent().map(|x| x as u64)
       }
    
       if *self == 0{
          return NTResult::Infinite;
       }
       let fctr = self.factor().unwrap();
       let mut result = 1;
       for i in 0..fctr.base.len(){
       if fctr.base[0] == 2 && fctr.power[0] > 2{
          let phi = 2u64<<(fctr.power[0]-1);
          result=phi;
       }
       else{
         let el = fctr.base[i];
         let phi =  (el.pow(fctr.power[i] as u32)/el)*(el-1);
         match result.lcm(phi){
             NTResult::Overflow => {return NTResult::Overflow;},
             NTResult::Eval(x) => result=x,
             _=> (),
         }
       }
       }
        return NTResult::Eval(result);
    }
    
    fn dedekind_psi(&self, k: Self) -> NTResult<Self> {
        if *self == 0{
         return NTResult::Infinite
       }
       
        let (k2, flag) = k.overflowing_shl(1);
        if flag {
            return NTResult::Overflow;
        }
        self.jordan_totient(k2).map(|y| y/self.jordan_totient(k).unwrap())
    }

   fn quadratic_residue(&self, n: Self) -> Self {
        if n == 0 {
            return self.wrapping_mul(*self)
        }
        ((*self as u128 * *self as u128) % n as u128) as Self
    }
    
    fn checked_quadratic_residue(&self, n: Self) -> NTResult<Self> {
        if n == 0 {
            return NTResult::from_option(self.checked_mul(*self),NTResult::Overflow)
        }
        NTResult::Eval(((*self as u128 * *self as u128) % n as u128) as Self)
    }
    

    fn product_residue(&self, other: Self, n: Self) -> Self {
        if n == 0 {
            return self.wrapping_mul(other)
        }
        ((*self as u128 * other as u128) % n as u128) as Self
    }
    
    fn checked_product_residue(&self, other: Self, n: Self) -> NTResult<Self> {
        if n == 0 {
            return NTResult::from_option(self.checked_mul(*self),NTResult::Overflow)
        }
        NTResult::Eval(((*self as u128 * other as u128) % n as u128) as Self)
    }

    fn exp_residue(&self, p: Self, modulus: Self) -> Self {
        
        if modulus == 0 {
            return self.pow(p as u32);
        }
        
        if modulus & 1 == 0 {
        
        let k = modulus.trailing_zeros() as u64;
        let s = modulus >> k;

        let reducer = (1 << k) - 1; // A shorthand for arithmetic over Z[2k]

        let k_rem = self.even_pow(p,reducer);

        let s_rem = self.odd_pow(p,s);
        
        let s_inv = s.inv_2()&reducer;
    
        let y = k_rem.wrapping_sub(s_rem).wrapping_mul(s_inv) & reducer;

        s_rem + s * y //)%modulus
    } else {
        self.odd_pow(p,modulus) //%modulus
    }
    }

    fn checked_exp_residue(&self, p: Self, modulus: Self) -> NTResult<Self> {
        if modulus == 0 {
            if  p > u32::MAX as u64 {
                return NTResult::Overflow;
            }
            match self.checked_pow(p as u32) {
                Some(x) => return NTResult::Eval(x),
                None => return NTResult::Overflow,
            };
        }

        NTResult::Eval(self.exp_residue(p,modulus))
    }

    fn legendre(&self, p: Self) -> i8 {
        let k = self.exp_residue((p - 1) >> 1, p);
        if k == 1 {
            return 1;
        };
        if k == p - 1 {
            return -1;
        };
        0i8
    }

    fn checked_legendre(&self, p: Self) -> NTResult<i8> {
        if p == 2 || !p.is_prime() {
            return NTResult::Undefined;
        }
        NTResult::Eval(self.legendre(p))
    }

    fn liouville(&self) -> i8 {
      if self.reducible(){
        return (*self as u32).liouville()
      }
        let primeomega = self.factor().unwrap().prime_omega();
        if primeomega & 1 == 0 {
            return 1;
        }
        -1
    }
    
    fn derivative(&self) -> NTResult<Self> {
       if *self < 94 {
         return (*self as u32).derivative().map(|y| y as Self)
       }
       
       let fctr = self.factor().unwrap();
       let mut sum : u64 = 0;
       
       for (power,base) in fctr.power_iter().zip(fctr.factor_iter()){
          match sum.checked_add((*power as u64)*(*self/base)){
            Some(x) => sum =x,
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
       0f64
    }
    
    fn mobius(&self) -> i8 {
     //if self.reducible(){
     //   return (*self as u32).mobius()
     // }
      let fctr = self.factor().unwrap();
     // if fctr.len() == 1{ // if only one factor then return -1
     //    return -1
     // }
      if !fctr.k_free(2){
         return 0;
      }
      let fctrsum = fctr.prime_omega();//fctr[1..].iter().step_by(2).sum::<Self>();
      if fctrsum&1 == 1{// if odd number of factors and square free
        return -1
      }
      1
    }

    fn jacobi(&self, k: Self) -> i8 {
        let mut n = *self;
        let mut p = k;
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

    fn checked_jacobi(&self, k: Self) -> NTResult<i8> {
        if k > 0 && k % 2 == 1 {
            return NTResult::Eval(self.jacobi(k));
        }
        NTResult::Undefined
    }
    
     fn kronecker(&self, k: Self) -> i8{
     let x = *self;
     if k == 0{
      if x == 1{
         return 1
      }
     return 0
    }
   if k == 1{
      return 1
   }
   let fctr = k.factor().unwrap();
   let mut start = 0;
   let mut res = 1;
   
   if fctr.base[0] ==  2{
     start = 1;
     if x&1 == 0{
     res = 0;
     }
     else if x % 8 == 1 || x % 8 == 7{
      res=1
     }
     else{
       res = (-1i8).pow(fctr.power[0] as u32)
     }
   }
   if fctr.base[0] == 2 && fctr.base.len() == 1{
     return res
   }
   for i in start..fctr.base.len(){
     res*=self.legendre(fctr.base[i]).pow(fctr.power[i] as u32);
   }
   res
}



  fn smooth(&self) -> NTResult<Self> {
       if *self == 0{
         return NTResult::Infinite
       }
       if *self == 1{
        return NTResult::DNE
       }
       self.factor().map(|x| x.max())
  
    }
    
    

    fn is_smooth(&self, b: Self) -> bool {
     match self.smooth(){
      NTResult::Infinite => false,
      NTResult::Eval(x) => x <= b, 
      _=> false,
     }
   }
   
    fn ord(&self, n: Self) -> NTResult<Self>{
      let ord_2 = |a: u64, p: u64| -> u64{
        let modulo = (1u64<<p)-1;
        let mut b = a&modulo;
   
         if b == 1{
            return 1;
         }
      for i in 1..p{
         b = b.wrapping_mul(b)&modulo;
         if b == 1{
           return 1<<i;
         }
    }
    return p;
    };

// Given ord(a,p)  calculate ord(a,p^n)
    let  pp_ord = |a: u64, b: u64, p: u64, e: u32| -> u64{
     for i in 0..e+1{
          if a.exp_residue(b*p.pow(i),p.pow(e)) ==1{
             return b*p.pow(i);
           }
   }
    return b*p.pow(e);
    };

    let p_ord = |a: u64, p: u64| -> u64{
   
   let fctr = (p-1).factor().unwrap();
   
   let mut m = p-1;
   for i in fctr.pair_iter(){
     for _ in 0..*i.1{
          if a.exp_residue(m/ *i.0,p) == 1{
            m = m/ *i.0;
          }
          else{
            break;
          }
     }
  }
 m
};


    if self.gcd(n) != 1{
       return NTResult::DNE;
    }
    let fctr = n.factor().unwrap();
    let mut fullord = 1u64;
    for i in fctr.pair_iter(){
     let mut ord : Self;
      if *i.0 == 2{
         ord = ord_2(*self,*i.1);
      }
      else{
        ord = p_ord(*self,*i.0);
        if *i.1 > 1{
           ord=pp_ord(*self,ord,*i.0,*i.1 as u32);
        }
      }
       fullord = fullord.lcm(ord).unwrap(); 
    }
    NTResult::Eval(fullord)
 }

}
