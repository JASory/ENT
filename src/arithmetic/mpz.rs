use crate::arithmetic::conversion::*;
use crate::arithmetic::inlineops::*;

use crate::arithmetic::sign::Sign;
use crate::arithmetic::sliceops::*;
use std::cmp::Ordering;
use crate::primitive::factorprim::poly_factor_128;
use crate::ntrait::NumberTheory;

/// Arbitrary precision integer (Z)
#[derive(Debug, Default, Clone, PartialEq)]
pub struct Mpz {
    pub(crate) sign: Sign,
    pub(crate) limbs: Vec<u64>,
}

impl std::fmt::Display for Mpz{
      fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
    
       let out = to_string(self.sign.clone(), self.limbs.clone());
       write!(f,"{}",out)
  }
}

impl std::convert::From<u64> for Mpz{
   fn from(val: u64) -> Self{
        Mpz::unchecked_new(Sign::Positive, vec![val])
   }
}

impl std::convert::From<u32> for Mpz{
   fn from(val: u32) -> Self{
        Mpz::unchecked_new(Sign::Positive, vec![val.into()])
   }
}

impl std::convert::From<i64> for Mpz{
   fn from(val: i64) -> Self{
        if val < 0i64 {
            return Mpz::unchecked_new(Sign::Negative, vec![val.unsigned_abs()]);
        }
        Mpz::unchecked_new(Sign::Positive, vec![val as u64])
      }  
}

impl std::convert::From<i32> for Mpz{
   fn from(val: i32) -> Self{
        if val < 0i32 {
            return Mpz::unchecked_new(Sign::Negative, vec![val.unsigned_abs() as u64]);
        }
        Mpz::unchecked_new(Sign::Positive, vec![val as u64])
      }  
}

impl std::convert::From<u128> for Mpz{
   fn from(val: u128) -> Self{
      let (x_lo, x_hi) = split(val);
        if x_hi == 0 {
            return Mpz::unchecked_new(Sign::Positive, vec![x_lo]);
        }
        Mpz::unchecked_new(Sign::Positive, vec![x_lo, x_hi])
   }
}

impl std::convert::From<i128> for Mpz{

   fn from(val: i128) -> Self {
        if val < 0i128 {
            let (x_lo, x_hi) = split(val.unsigned_abs());
            if x_hi == 0 {
                return Mpz::unchecked_new(Sign::Negative, vec![x_lo]);
            }
            Mpz::unchecked_new(Sign::Negative, vec![x_lo, x_hi])
        } else {
            Mpz::from(val as u128)
        }
    }
  }

impl std::convert::TryFrom<Mpz> for u64{
      type Error = &'static str;

   fn try_from(x: Mpz) ->  Result<Self, Self::Error>{
        
        match x.len() {
         0 => Ok(0u64),
         1 => Ok(x.limbs[0]),
         _=> Err("Multiprecision value exceeds 2^64"),
        }
   }

}  

impl std::convert::TryFrom<Mpz> for i64{
     type Error = &'static str;

  fn try_from(x: Mpz) -> Result<Self,Self::Error>{
     
      match x.len(){
          0 => Ok(0i64),
          1 => {
             let value = x.limbs[0];
            if (value>>63)==1{
               return  Err("Single precision value exceeds 2^63");
             }
            let mut res : i64 = value as i64;
            if x.sign == Sign::Negative{
               res = -res;
            }
            return Ok(res);
           }
          _=> Err("Multiprecision value exceeds 2^64"),
      }
  }
}  
  /*
 impl std::str::FromStr for Mpz{
 
 } 

*/
impl Mpz {
    /**

    ```
     use number_theory::Mpz;  // includes arbitrary-precision arithmetic
     use number_theory::Sign; // allows sign
     use number_theory::NumberTheory; // includes methods from NumberTheory trait

     let bignum = Mpz::new(Sign::Positive, vec![5,5]);

     let fortyfour = Mpz::from(44u128);

     let fortyfour_neg = Mpz::from(-44i128);

     let twopow = Mpz::from(128u64);


     assert_eq!("92233720368547758085", bignum.to_string());
     assert_eq!("44", fortyfour.to_string());
     assert_eq!("-44", fortyfour_neg.to_string());
     assert_eq!("128", twopow.to_string());

    ```

    */
      /// Returns a positive x from a big-endian vector representation 
    pub fn u_new(limbs: Vec<u64>) -> Self {
        Mpz::from_slice(Sign::Positive, &limbs[..])
    }
      /// Returns a x from a big-endian vector representation and a sign 
    pub fn new(sign: Sign, limbs: Vec<u64>) -> Self {
        Mpz::from_slice(sign, &limbs[..])
    }
     /// Returns a positive x from a big-endian vector, this function does not eliminate extraneous leading zeros 
    pub fn unchecked_u_new(limbs: Vec<u64>) -> Self {
        Mpz {
            sign: Sign::Positive,
            limbs,
        }
    }

  /// Returns a  x from a big-endian vector and sign, this function does not eliminate extraneous leading zeros 
    pub fn unchecked_new(sign: Sign, limbs: Vec<u64>) -> Self {
        Mpz { sign, limbs }
    }
/// Returns x from a big_endian slice and sign
    pub fn from_slice(sign: Sign, x: &[u64]) -> Self {
        let mut limbs = vec![];
        limbs.extend_from_slice(&x[..sig_pos(x)]);
        if limbs.is_empty() {
            limbs.push(0u64)
        }
        Mpz { sign, limbs }
    }
    
    /// Converts to 32-bit unsigned integer if it can fit 
    pub fn to_u32(&self) ->  Option<u32>{
       match self.len() {
            0 => Some(0u32),
            1 => {
                   if self.limbs[0] > u32::MAX as u64 {
                      return None
                   } 
                   Some(self.limbs[0] as u32)
                 }
            _ => None,
        }
    }
    
     /// Converts to 128-bit unsigned integer if it can fit
    pub fn to_u128(&self) -> Option<u128> {
        match self.len() {
            0 => Some(0u128),
            1 => Some(self.limbs[0] as u128),
            2 => Some(fuse(self.limbs[1], self.limbs[0])),
            _ => None,
        }
    }
    /// Returns the vector representation of n
    pub fn to_vector(&self) -> Vec<u64> {
        self.limbs.clone()
    }
   
     /// Set an integer between 2^(n-1);2^n
    pub fn rand(len: usize) -> Self {
      let mut k = len/64;
      if len%64 != 0{
        k+=1;
      }
      let mut interim = (0..k).map(|_| u64::rng()).collect::<Vec<u64>>();
      
      interim[k-1] &= (1<<(len%64))-1;
      let mut result = Mpz::unchecked_new(Sign::Positive, interim);
      result.set_bit(len-1);
      result
      
      
    }
    
    
    /// Creates a integer from limbs generated by a function 
    pub fn set_limbs(len: usize, gen: fn() -> u64) -> Self {
        let interim = (0..len).map(|_| gen()).collect::<Vec<u64>>();
        Mpz::unchecked_new(Sign::Positive, interim)
    }
    /// Converts n to a radix-10 string
    pub fn to_string(&self) -> String {
        to_string(self.sign.clone(), self.limbs.clone())
    }
    /// Converts n to a radix-16 string in linear time
    pub fn to_hex_string(&self) -> String{
       to_hex_string(self.sign.clone(),self.limbs.clone())
    }
    // temporary placeholder function to flip negative zeros
    pub(crate) fn fix_zero(&mut self) {
        if self.len() == 0 {
            self.limbs.push(0u64)
        }
        if self.len() == 1 && self.limbs[0] == 0 && self.sign == Sign::Negative {
            self.sign = Sign::Positive;
        }
    }
    /**
      Returns the polynomial representation of self in the form of the coefficient vector. Here we see that 50620 to radix-127 is 3x^2 + 17x + 74.
      Accepts all radix in the interval 0;2^64-1.

    ```
    use number_theory::Mpz;

     let value = Mpz::from(50620u64);
     assert_eq!(value.to_radix_vec(127),vec![74,17,3])
    ```
    */

    pub fn to_radix_vec(&self, radix: u64) -> Vec<u64> {
        let mut k = vec![];
        let mut x = self.clone().limbs;
        loop {
            let idx = sig_pos(&x[..]);

            if idx == 0 {
                break;
            }

            k.push(div_slice(&mut x[..idx], radix));
        }
        k.push(x[0] % radix);
        remove_lead_zeros(&mut k);
        k
    }

    /**
      Conversion from radix-10^n string.

    ```
    use number_theory::Mpz;
     let num = "-3141592653589793238462643383279502884197169399375105820974944592307816406286".to_string();
     let machine = Mpz::from_string(&num).unwrap();

     assert_eq!(machine.to_string(),num)
    ```
    */

    pub fn from_string(x: &str) -> Option<Self> {
        if x.is_empty() {
            return None;
        }

        let ch = x.chars().next().unwrap();
        let mut sign = Sign::Positive;
        let mut k = x;
        if ch == '-' {
            sign = Sign::Negative;
            let mut chars = x.chars();
            chars.next();
            k = chars.as_str();
        }

        if ch == '+' {
            let mut chars = x.chars();
            chars.next();
            k = chars.as_str();
        }

        from_string(k).map(|y| Mpz::unchecked_new(sign, y))
    }
    
    /// Returns 0
    pub fn zero() -> Self {
        Mpz::unchecked_new(Sign::Positive, vec![0])
    }
    /// Returns positive 1
    pub fn one() -> Self {
        Mpz::unchecked_new(Sign::Positive, vec![1])
    }
    /// Returns positive 2
    pub fn two() -> Self{
       Mpz::unchecked_new(Sign::Positive, vec![2])
    }
    /// Additive inverse of n, evaluated in-place 
    pub fn neg(&mut self) {
        self.sign = self.sign.neg();
    }
     /// Returns the absolute value of n  
    pub fn abs(&self) -> Self {
        Self::new(Sign::Positive, self.limbs.clone())
    }
    
    /// Checks if n == 0
    pub fn is_zero(&self) -> bool {
       if self.len() == 0{
         return true
       }
       if self.len() == 1 && self.limbs[0] == 0{
         return true
       }
       false
    }
    /// Checks if n >= 0
    pub fn is_positive(&self) -> bool {
        self.sign == Sign::Positive
    }
     /// Checks if n in 2Z
    pub fn is_even(&self) -> bool {
        self.limbs[0] & 1 == 0
    }
    /// Checks if n in 2^Z
    pub fn is_power_of_two(&self) -> bool {
        if self.len() < 2 {
            return self.to_u128().unwrap().is_power_of_two();
        }
        if !self.lead_digit().is_power_of_two() {
            return false;
        }
        for i in self.limbs[..self.len() - 1].iter() {
            if *i != 0 {
                return false;
            }
        }
        true
    }

    pub(crate) fn is_fermat(&self) -> bool {
        if self.limbs[0] != 1 {
            return false;
        }

        let lead = self.limbs[..].last().unwrap();

        let mut flag = 0u64;
        for i in 0..64 {
            if *lead == 1u64 << i {
                flag = 1; // set flag
                break; // end loop
            }
        }

        if flag == 0 {
            return false;
        } // if the flag is not set return false

        for i in self.limbs[1..self.len() - 1].iter() {
            if *i != 0u64 {
                return false;
            }
        }
        true
    }

    pub(crate) fn is_mersenne(&self) -> Option<u64> {
        let mut flag = 0u64;
        let mut start = 1u64;
        let lead = self.limbs.last().unwrap();
        for i in 0..64 {
            if *lead == start {
                flag = i;
                break;
            }
            start = (start << 1) + 1;
        }
        if flag == 0 {
            return None;
        }
        for i in self.limbs[..self.len() - 2].iter() {
            if *i != u64::MAX {
                return None;
            }
        }

        Some(flag + 1 + 64u64 * (self.len() - 1) as u64)
    }
    /// Sets the bit at 2^k to 1, if it is already 1 then this does not change
    pub fn set_bit(&mut self, k: usize) {
        self.limbs[k / 64usize] |= 1 << (k % 64)
    }
    /// Flips the bit at 2^k
    pub fn flip_bit(&mut self, k: usize){
         self.limbs[k / 64usize] ^= 1 << (k % 64)
    }
   /// Check if the bit in the k-place is set  
    pub fn check_bit(&self, index: usize) -> bool {
        self.limbs[index / 64usize] & (1 << (index % 64)) > 0
    }
   /// Change the sign   
    pub fn set_sign(&mut self, sign: Sign) {
        self.sign = sign
    }
    /// Returns the number of 64-bit machine words used in this representation
    pub fn len(&self) -> usize {
        self.limbs.len()
    }
    /// Returns the lead digit in 2^64 representation
    pub fn lead_digit(&self) -> u64 {
        *self.limbs[..].last().unwrap()
    }
   /// Returns k such that  d*2^k = x
    pub fn trailing_zeros(&self) -> u64 {
        // Trailing zeros
        let mut idx: u64 = 0;

        for i in self.limbs.iter() {
            if i == &0u64 {
                idx += 1;
            } else {
                break;
            }
        }
        if idx == self.len() as u64 {
            64u64 * idx
        } else {
            self.limbs[idx as usize].trailing_zeros() as u64 + 64u64 * idx
        }
    }
    /// Returns the number of extraneous zeros in bit representation
    pub fn leading_zeros(&self) -> u64 {
        let mut idx: u64 = 0;

        for i in self.limbs.iter().rev() {
            if i == &0u64 {
                idx += 1;
            } else {
                break;
            }
        }
        if idx == self.len() as u64 {
            64u64 * idx
        } else {
            self.limbs[self.len() - 1 - idx as usize].leading_zeros() as u64 + 64u64 * idx
        }
    }
  /// Returns the number of bits used to represent x 
    pub fn bit_length(&self) -> u64 {
        64u64 * self.len() as u64 - self.leading_zeros()
    }
   /// Removes extraneous zeros 
    pub fn normalize(&mut self) {
        remove_lead_zeros(&mut self.limbs);
        if self.len() == 0 {
            self.limbs.push(0u64)
        };
    }

    /*
    Equality Operations

    */
  /// Compares the magnitude of x and y, ignoring sign
    pub fn u_cmp(&self, other: &Self) -> Ordering {
        cmp_slice(&self.limbs[..], &other.limbs[..])
    }
   /// Checks if x mod n === c 
    pub fn congruence_u64(&self, n: u64, c: u64) -> bool {
        let mut interim = mod_slice(&self.limbs[..], n);
        if self.sign == Sign::Negative {
            interim = n - interim;
        }
        interim == c
    }

    pub(crate) fn add_modinv(&self, n: &Self) -> Self {
        // additive modular inverse
        let mut k = n.clone();
        let mut selfie = self.clone();

        if self.u_cmp(n) == Ordering::Greater {
            selfie = selfie.ref_euclidean(n).1;
        }

        sub_slice(&mut k.limbs, &selfie.limbs);
        k.normalize();
        k.sign = Sign::Positive;
        k
    }

    /**

       ```
             use number_theory::Mpz;
         let mut one = Mpz::one();

         one.successor();  //Applies the successor function ignoring sign

         assert_eq!("2", one.to_string())

       ```


    */

    pub fn successor(&mut self) {
        if self.len() == 0 {
            self.limbs.push(1)
        }
        let mut carry = 1u8;
        for i in self.limbs.iter_mut() {
            carry = adc(carry, *i, 0, i);
            if carry == 0 {
                break;
            }
        }
        if carry > 0u8 {
            self.limbs.push(1u64)
        }
    }
     /// Increments by n amount
    pub fn inc_by(&mut self, n: u64){
      if self.len() == 0 {
            self.limbs.push(n)
        }
        
        add_slice(&mut self.limbs[..],&[n]);
    }
    
    /// Inverse successor function. Subtracts 1 ignoring sign
    pub fn predecessor(&mut self) {
        if self.len() == 0 {
            panic!("Value already zero")
        }
        let mut carry = 1u8;
        for i in self.limbs.iter_mut() {
            carry = sbb(carry, *i, 0, i);
            if carry == 0 {
                break;
            }
        }
        if carry > 0u8 {
            panic!("Value already zero")
        }
    }

    /*
      Precursor NT functions, unsigned
    */

    pub(crate) fn u_quadratic_residue(&self, n: &Self) -> Self {
        self.ref_product(self).ref_euclidean(n).1
    }

    pub(crate) fn u_mul_mod(&self, other: &Self, n: &Self) -> Self {
        self.ref_product(other).ref_euclidean(n).1
    }

    pub(crate) fn u_mod_pow(&self, y: &Self, modulo: &Self) -> Self {
        if modulo == &Mpz::zero() {
            return Mpz::zero();
        }

        let mut z = Mpz::one();
        let mut base = self.clone().ref_euclidean(modulo).1;
        let one = Mpz::one();

        let mut pow = y.clone();

        if pow == Mpz::zero(){
            return z;
        }

        if pow == Mpz::one() {
            return base;
        }

        while pow.u_cmp(&one) == Ordering::Greater {
            if pow.len() == 0 {
                break;
            }

            if pow.is_even() {
                base = base.u_quadratic_residue(modulo);
                pow.mut_shr(1);
                remove_lead_zeros(&mut pow.limbs);
            } else {
                z = base.u_mul_mod(&z, modulo);
                remove_lead_zeros(&mut base.limbs);
                base = base.u_quadratic_residue(modulo);

                sub_slice(&mut pow.limbs[..], &one.limbs[..]);
                pow.mut_shr(1);
                remove_lead_zeros(&mut pow.limbs);
            }
        }
         base.u_mul_mod(&z, modulo)
    }

    /// Returns self^x
    pub fn pow(&self, x: u64) -> Self {
        let mut z = Mpz::one();
        let mut base = self.clone();
        let mut pow = x;

        if pow == 0 {
            return z;
        }

        while pow > 1 {
            if pow % 2 == 0 {
                base = base.ref_product(&base);
                pow >>= 1;
            } else {
                z = base.ref_product(&z);
                base = base.ref_product(&base);
                pow = (pow - 1) >> 1;
            }
        }

        base.ref_product(&z)
    }

    pub(crate) fn word_div(&self, x: u64) -> (Self, u64) {
        let mut quotient = self.clone();
        let remainder = div_slice(&mut quotient.limbs[..], x);
        (quotient, remainder)
    }

    fn delta(&self, other: &Self) -> Self {
        if self.u_cmp(other) == Ordering::Greater {
            let mut k = self.clone();
            sub_slice(&mut k.limbs[..], &other.limbs[..]);
            return k;
        }
        if self.u_cmp(other) == Ordering::Less {
            let mut k = other.clone();
            sub_slice(&mut k.limbs[..], &self.limbs[..]);
            k
        } else {
            Mpz::zero()
        }
    }

    fn mod_sqr_1(&self, n: &Self) -> Self {
        let mut k = self.ref_product(self);
        k.successor();
        k.ref_euclidean(n).1
    }

    fn gcd(&self, other: &Self) -> Self {
        let mut a = self.clone();
        let mut b = other.clone();

        while b != Mpz::zero() {
            let t = b.clone();

            b = a.ref_euclidean(&b).1;

            b.normalize();
            a = t;
        }
        a
    }
    
// FIXME 
    pub(crate) fn rho_mpz(&self) -> Self {
    
     match self.to_u128(){
      Some(x) => Mpz::from(poly_factor_128(x)),
      None => {
        let mut x = Mpz::two();
        let mut y = Mpz::two();
        let mut d = Mpz::one();
        loop {
            while d.is_unit(){
                x = x.mod_sqr_1(self);
                y = y.mod_sqr_1(self).mod_sqr_1(self).ref_euclidean(self).1;
                d = x.delta(&y).gcd(self);
            }

            if d.is_prime() {
                return d;
            }
            x = Mpz::from(u64::rng());
            y = x.clone();
            d = Mpz::one();
        }
       
       },
    }
    }

    /// Approximation of the natural logarithm ln(x)
    pub fn ln(&self) -> f64 {
        let mut iter_sqrt = 1u64;
        let mut n = self.sqrt().0;
        while n.len() > 1 {
            n = n.sqrt().0;
            iter_sqrt += 1;
        }

        (u64::try_from(n).unwrap() as f64).ln() * 2f64.powf(iter_sqrt as f64)
    }
    /// Approximation of the base-2 logarithm
    pub fn log2(&self) -> f64 {
        self.ln() * std::f64::consts::LOG2_E
    }
    /// Approximation of the base-10 logarithm
    pub fn log10(&self) -> f64 {
        self.ln() * 0.43429448190325176
    }
    /// Approximation of the base-k logarithm, where k is a Real.
    pub fn log(&self, log: f64) -> f64 {
        self.ln() * log.ln().recip()
    }
    /// Iterated logarithm
    pub fn iter_log(&self, log: f64) -> u8 {
        let mut first_log = self.log(log);
        let mut count = 1u8;
        // 1.444667861009766

        while first_log > 1.0 {
            first_log = first_log.log(log);
            count += 1
        }
        count
    }
    
    /// Finite ring variant of Extended Euclidean algorithm, see trait definition of extended gcd
    pub fn eea(&self, other: Self) -> (Self, Self, Self) {
    
    let (gcd, bezout_1, bezout_2) = self.extended_gcd(other.clone());
    
    (gcd, bezout_1.residue(other.clone()), bezout_2.residue(self.clone()))
     
    }

    /**
    Sequential Interval Residue Product, a generalization of the k-factorials
    
      ```

         // see the sirp crate for greater extension of this functionality
          use number_theory::Mpz;
          
         let factorial_100 = Mpz::sirp(1,100,1,0);
         
         // Even. Odd double factorial would be sirp(1,100,2,1)
         let doublefact_100 = Mpz::sirp(1,100,2,0);
         assert_eq!(
         "3424322470251197624824643289520818597\
          5118675053719198827915654463488000000000000",doublefact_100.to_string())
      ```
    */

    pub fn sirp(infimum: u64, supremum: u64, modulo: u64, residue: u64) -> Self {
        let mut sirp = Mpz::one();
        let mut acc = 1u64; // accumulator for factors

        for i in infimum..supremum + 1 {
            // inclusive range

            if i % modulo == residue {
                // if i is of the residue class modulo n then multiply   If you set residue as stop mod n then you get the k-factorials

                if i >= 4294967296 {
                    acc = i
                }
                if acc < 4294967296 {
                    acc *= i
                }
                if acc >= 4294967296 {
                    let carry = scale_slice(&mut sirp.limbs[..], acc);

                    if carry > 0 {
                        sirp.limbs.push(carry)
                    }
                    acc = 1u64;
                } // end if
            } // else
        }

        let carry = scale_slice(&mut sirp.limbs[..], acc);
        if carry > 0 {
            sirp.limbs.push(carry)
        }
        sirp
    }
    /**
        Conditional Interval Product computes the product of integers satisfying an unary function. In this example the unary function is the
         primality function. In practice it can be any function that borrows a u64 and returns a boolean. 

      ```

          use number_theory::Mpz;
          use number_theory::NumberTheory;

         let primorial_25 = Mpz::cip(1,100,u64::is_prime);


      ```
    */

    // conditional interval product
    pub fn cip(infimum: u64, supremum: u64, cond: fn(&u64) -> bool) -> Self {
        let mut cip = Mpz::one();
        let mut acc = 1u64; // accumulator for factors
        for i in infimum..supremum + 1 {
            if cond(&i) {
                if i >= 4294967296 {
                    acc = i
                }
                if acc < 4294967296 {
                    acc *= i
                }
                if acc >= 4294967296 {
                    let carry = scale_slice(&mut cip.limbs[..], acc);

                    if carry > 0 {
                        cip.limbs.push(carry)
                    }
                    acc = 1u64;
                } // end if
            } // else
        }

        let carry = scale_slice(&mut cip.limbs[..], acc);
        if carry > 0 {
            cip.limbs.push(carry)
        }
        cip
    }
}
