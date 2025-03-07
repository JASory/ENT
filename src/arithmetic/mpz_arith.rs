use crate::arithmetic::inlineops::*;
use crate::arithmetic::muldiv::*;
use crate::arithmetic::sliceops::*;
use crate::Mpz;
use crate::arithmetic::sign::Sign;
use std::cmp::Ordering;
use crate::ntrait::NumberTheory;

impl std::cmp::PartialOrd for Mpz{
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        match (self.sign,other.sign){
         (Sign::Negative,Sign::Positive) => Some(std::cmp::Ordering::Less),
         (Sign::Positive,Sign::Negative) => Some(std::cmp::Ordering::Greater),
         _=> Some(self.u_cmp(other)),
       }
     }  
}

impl std::cmp::Eq for Mpz{}

impl std::cmp::Ord for Mpz{
    fn cmp(&self,other: &Self) -> std::cmp::Ordering{
       match (self.sign,other.sign){
         (Sign::Negative,Sign::Positive) => std::cmp::Ordering::Less,
         (Sign::Positive,Sign::Negative) => std::cmp::Ordering::Greater,
         _=> self.u_cmp(other),
       }
        
    }
}

impl Mpz {
    /*
      Shifting operations

    */
   /// Performs bitwise AND operation between x and y storing the result in x
    pub fn mut_and(&mut self, other: &Self) {
        for (i, j) in self.limbs.iter_mut().zip(other.limbs.iter()) {
            *i &= j
        }
        self.limbs.truncate(other.len());
        self.normalize();
    }

   /// Performs bitwise OR operation between x and y storing the result in x
    pub fn mut_or(&mut self, other: &Self) {
        for (i, j) in self.limbs.iter_mut().zip(other.limbs.iter()) {
            *i |= j
        }

        if self.len() < other.len() {
            self.limbs.extend_from_slice(&other.limbs[self.len()..])
        }
    }
   /// Performs bitwise XOR operation between x and y storing the result in x
    pub fn mut_xor(&mut self, other: &Self) {
        for (i, j) in self.limbs.iter_mut().zip(other.limbs.iter()) {
            *i ^= j
        }

        if self.len() < other.len() {
            self.limbs.extend_from_slice(&other.limbs[self.len()..])
        }
    }
  /// Performs bitwise NOT on x in-place
    pub fn mut_not(&mut self) {
        for i in self.limbs.iter_mut() {
            *i = !*i
        }
        self.normalize();
    }
/// Shift-left k places in-place
    pub fn mut_shl(&mut self, shift: usize) {
        //let mut k = self.clone();
        let mut trail_zeroes = vec![0; shift / 64usize];

        let carry = shl_slice(&mut self.limbs[..], (shift % 64usize) as u32);

        trail_zeroes.extend_from_slice(&self.limbs[..]);

        if carry > 0 {
            trail_zeroes.push(carry)
        }

        self.limbs = trail_zeroes;
    }
    
/// Shift-left k places in-place
    pub fn mut_shr(&mut self, shift: usize) {
        let mut carry = 0u64;

        let mut vector: Vec<u64> = self
            .limbs
            .drain((shift / 64usize)..self.limbs.len())
            .collect();
        let sub_shift = shift % 64usize;

        for i in vector.iter_mut().rev() {
            carry = carry_shr(carry, *i, sub_shift as u32, i);
        }

        self.limbs = vector;
    }

/// Returns x * 2^k
    pub fn shl(&self, shift: usize) -> Mpz {
        let mut k = self.clone();
        k.mut_shl(shift);
        k
    }
/// Returns x / 2^k
    pub fn shr(&self, shift: usize) -> Mpz {
        let mut k = self.clone();
        k.mut_shr(shift);
        k
    }
/// Returns x AND y 
    pub fn and(&self, other: &Self) -> Self {
        let mut k = self.clone();
        k.mut_and(other);
        k
    }
///  Returns x OR y 
    pub fn or(&self, other: &Self) -> Self {
        let mut k = self.clone();
        k.mut_or(other);
        k
    }
/// Returns x XOR y
    pub fn xor(&self, other: &Self) -> Self {
        let mut k = self.clone();
        k.mut_xor(other);
        k
    }

    /*
        Arithmetic operations
    */
   /// x+y stored in x
    pub fn mut_addition(&mut self, mut other: Self) {
        let carry: u8;

        if self.sign == other.sign {
            if self.limbs.len() < other.limbs.len() {
                let len = self.len();
                self.limbs.extend_from_slice(&other.limbs[len..]);
                carry = add_slice(&mut self.limbs[..], &other.limbs[..len]);
            } else {
                carry = add_slice(&mut self.limbs[..], &other.limbs[..]);
            }
            if carry == 1u8 {
                self.limbs.push(1u64)
            }
        } else {
            if self.u_cmp(&other) == Ordering::Less {
                sub_slice(&mut other.limbs[..], &self.limbs[..]);
                *self = other;
                self.normalize();
            } else if self.u_cmp(&other) == Ordering::Equal {
                self.limbs.truncate(0);
                self.limbs.push(0);
                self.sign = Sign::Positive;
            } else {
                sub_slice(&mut self.limbs[..], &other.limbs[..]);
                self.normalize();
            }
        }
    }
   /// Returns x+y 
    pub fn ref_addition(&self, other: &Self) -> Self {
        let mut k = self.clone();
        k.mut_addition(other.clone());
        k
    }
   /// x-y stored in x 
    pub fn mut_subtraction(&mut self, mut other: Self) {
        other.neg();
        self.mut_addition(other)
    }
    /// Returns x-y
    pub fn ref_subtraction(&self, other: &Self) -> Self {
        let mut k = self.clone();
        k.mut_subtraction(other.clone());
        k
    }
   /// Returns x*y 
    pub fn ref_product(&self, other: &Self) -> Self {
        if self == &Mpz::zero() {
            return Mpz::zero();
        }

        if other == &Mpz::zero() {
            return Mpz::zero();
        }

        if self.is_unit() {
            return Mpz::unchecked_new(self.sign.mul(&other.sign), other.limbs.clone());
        }

        if other.is_unit() {
            return Mpz::unchecked_new(self.sign.mul(&other.sign), self.limbs.clone());
        }

        let mut t = vec![0u64; self.len() + other.len() + 1];

        mul_slice(&self.limbs[..], &other.limbs[..], &mut t[..]);
        remove_lead_zeros(&mut t);
        Mpz::unchecked_new(self.sign.mul(&other.sign), t)
    }

/// x*y stored in x 
    pub fn mut_product(&mut self, other: Self) {
        if self == &Mpz::zero() {}

        if other == Mpz::zero() {
            self.limbs.truncate(0);
            self.limbs.push(0u64);
            self.sign = Sign::Positive;
        }

        if self.is_unit() {
            self.limbs.truncate(0);
            self.limbs.extend_from_slice(&other.limbs[..]);
            self.sign = self.sign.mul(&other.sign);
        }

        if other.is_unit() {
            self.sign = self.sign.mul(&other.sign);
        } else {
            let mut t = vec![0u64; self.len() + other.len() + 1];

            mul_slice(&self.limbs[..], &other.limbs[..], &mut t[..]);
            remove_lead_zeros(&mut t);
            self.limbs = t;
            self.sign = self.sign.mul(&other.sign);
        }
    }
    /// Unsigned in-place multiply by x and add y 
    pub fn mut_scale_add(&mut self, x: u64, y: u64){
     let mut carry = scale_slice(&mut self.limbs[..],x);
     carry += add_slice(&mut self.limbs[..],&[y]) as u64;
     if carry > 0{
       self.limbs.push(carry);
     }
    }
    
     /// Unsigned in-place multiply by x and add y 
    pub fn scale_add(&self, x: u64, y: u64) -> Self{
     let mut k = self.clone();
     k.mut_scale_add(x,y);
     k
    }
    
    
    /// x*x squaring 
    pub fn sqr(&self) -> Self {
      self.ref_product(self)
    }
    
    /// x*x squaring in-place
    pub fn mut_sqr(&mut self){
      self.mut_product(self.clone());
    }
        
  /// Returns x/y and x mod y 
    pub fn ref_euclidean(&self, other: &Self) -> (Self, Self) {
        let mut dividend = Mpz::from_slice(Sign::Positive, &self.limbs[..]);

        if dividend == Mpz::zero() {
            return (Mpz::zero(), Mpz::zero());
        }

        if dividend.len() == 0usize {
            return (Mpz::zero(), Mpz::zero());
        }

        if other.len() == 0usize || other == &Mpz::zero() {
            panic!("Division by zero is undefined in Z")
        }

        if other.len() == 1 {
            if other.limbs == [1] {
                return (dividend, Mpz::zero());
            }

            let rem = div_slice(&mut dividend.limbs, other.limbs[0]);

            remove_lead_zeros(&mut dividend.limbs);
            return (dividend, Mpz::unchecked_new(Sign::Positive, vec![rem]));
        }

        if dividend.u_cmp(other) == Ordering::Equal {
            return (Mpz::one(), Mpz::zero());
        }

        if dividend.u_cmp(other) == Ordering::Less {
            return (Mpz::zero(), dividend);
        }

        let shift = other.limbs.last().unwrap().leading_zeros() as usize;

        if shift == 0 {
            let (quo, rem) = euclidean_slice(&mut dividend.limbs, &other.limbs[..]);
            (
                Mpz::unchecked_new(Sign::Positive, quo),
                Mpz::unchecked_new(Sign::Positive, rem),
            )
        } else {
            let (q, r) =
                euclidean_slice(&mut dividend.shl(shift).limbs, &other.shl(shift).limbs[..]);

            (
                Mpz::unchecked_new(Sign::Positive, q),
                Mpz::unchecked_new(Sign::Positive, r).shr(shift),
            )
        }
    }
}
