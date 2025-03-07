use crate::ntrait::MontArith;

impl MontArith for u32{

   fn to_mont(self, ring: Self) -> Self{
      (((self as u64)<<32) % (ring as u64)) as u32
   }
   
   fn inv_2(self) -> Self{
     let mut est = machine_prime::INV_8[((self>>1) & 0x7F) as usize] as u32;
     est = 2u32.wrapping_sub(est.wrapping_mul(self)).wrapping_mul(est);
     est = 2u32.wrapping_sub(est.wrapping_mul(self)).wrapping_mul(est);
     est
   }
   
   fn n_identity(self) -> Self{
      (u32::MAX%self).wrapping_add(1)
   }
   
   fn mont_sub(self, subtrahend: Self, ring: Self) -> Self{
       if subtrahend > self {
           return ring.wrapping_sub(subtrahend.wrapping_sub(self));
       }
      self.wrapping_sub(subtrahend)
   }
   
   fn to_z(self,inv: Self, ring: Self) -> Self{
      let lo = self.wrapping_mul(inv);
      let lo = ((lo as u64).wrapping_mul(ring as u64)>>32) as u32;

      lo.wrapping_neg().wrapping_add(ring)
    
   }
   
   fn mont_add(self, addend: Self, ring: Self) -> Self{
       unimplemented!()
   }
   
   fn mont_prod(self, multiplicand: Self, inv: Self, ring: Self) -> Self{
      let interim = (self as u64).wrapping_mul(multiplicand as u64);
      let lo = (interim as u32).wrapping_mul(inv);
      let lo = ((lo as u64).wrapping_mul(ring as u64)>>32) as u32;
      let hi = (interim>>32)  as u32;
    
      if hi < lo{
         hi.wrapping_sub(lo).wrapping_add(ring)
      }
      else{
         hi.wrapping_sub(lo)
      }

   }
  
   fn mont_sqr(self, inv: Self, ring: Self) -> Self{
      self.mont_prod(self,inv,ring)
   }
  
   fn mont_pow(self,mut one: Self,mut p: Self, inv: Self, ring: Self) -> Self{
     let mut base = self;
     while p > 1 {
        if p & 1 == 0 {
            base = base.mont_sqr(inv, ring);
            p >>= 1;
        } else {
            one = one.mont_prod(base, inv, ring);
            base = base.mont_sqr(inv,ring);
            p = (p - 1) >> 1;
        }
    }
    one.mont_prod(base,inv,ring)
   }
   
    fn odd_pow(self, p: Self, n: Self) -> Self{
      
    let one = n.n_identity();
    let base = self.to_mont(n);
    let npi = n.inv_2();
    base.mont_pow(one,p,npi,n).to_z(npi,n)
   }
   
   fn even_pow(self, p: Self,n: Self) -> Self{
       let mut z = 1;
       let mut base = self;
       let mut pow = p;
       while pow > 1 {
           if pow & 1 == 0 {
               base = base.wrapping_mul(base)&n;
                pow >>= 1;
            } else {
                 z = base.wrapping_mul(z)&n;
              base = base.wrapping_mul(base)&n;
               pow = (pow - 1) >> 1
        }
    }
      base.wrapping_mul(z)&n
      } 
   
}


impl MontArith for u64{

   fn to_mont(self, ring: Self) -> Self{
        machine_prime::to_mont(self,ring)
   }
   
   fn inv_2(self) -> Self{
     machine_prime::mul_inv2(self)
   }
   
   fn n_identity(self) -> Self{
      machine_prime::one_mont(self)   
   }
   
   fn mont_sub(self, subtrahend: Self, ring: Self) -> Self{
     machine_prime::mont_sub(self,subtrahend,ring)
   }
   
   fn mont_add(self, addend: Self, ring: Self) -> Self{
      unimplemented!()
   }
   
   fn to_z(self,inv: Self, ring: Self) -> Self{
      let lo = self.wrapping_mul(inv);
      let lo = ((lo as u128).wrapping_mul(ring as u128)>>64) as u64;

      lo.wrapping_neg().wrapping_add(ring)
    
   }
   
   fn mont_prod(self, multiplicand: Self, inv: Self, ring: Self) -> Self{
      machine_prime::mont_prod(self,multiplicand,ring,inv)
   }
  
   fn mont_sqr(self, inv: Self, ring: Self) -> Self{
      self.mont_prod(self,inv,ring)
   }
  
   fn mont_pow(self,one: Self,p: Self, inv: Self, ring: Self) -> Self{
       machine_prime::mont_pow(self,one,p,ring,inv)
   }
   
    fn odd_pow(self, p: Self, n: Self) -> Self{
      
    let one = n.n_identity();
    let base = self.to_mont(n);
    let npi = n.inv_2();
    if base < 2{
       return base;
    }
    base.mont_pow(one,p,npi,n).to_z(npi,n)
   }
   
   fn even_pow(self, mut pow: Self,n: Self) -> Self{
       let mut z = 1;
       let mut base = self;
       while pow > 1 {
           if pow & 1 == 0 {
               base = base.wrapping_mul(base)&n;
                pow >>= 1;
            } else {
                 z = base.wrapping_mul(z)&n;
              base = base.wrapping_mul(base)&n;
               pow = (pow - 1) >> 1
        }
    }
      base.wrapping_mul(z)&n
      }

 }  

#[inline(always)]
const fn split_to_u128(x: u128) -> (u128, u128) { 
    let (lo,hi) = unsafe { std::mem::transmute::<u128, (u64, u64)>(x) };
    (hi as u128, lo as u128)
}



impl MontArith for u128{

   fn to_mont(self, ring: Self) -> Self{
      machine_prime::to_mont_128(self%ring,ring)
  }
   
   fn inv_2(self) -> Self{
       machine_prime::mul_inv2_128(self)
   }
   
   fn n_identity(self) -> Self{
      machine_prime::one_mont_128(self)   
   }
   
   fn to_z(self, inv: Self, ring: Self) -> Self{
         let lo = self.wrapping_mul(inv);
         let lo = machine_prime::u256prod_lo(lo,ring);

         lo.wrapping_neg().wrapping_add(ring)
   }
   
   
   fn mont_sub(self, subtrahend: Self, ring: Self) -> Self{
            machine_prime::mont_sub_128(self,subtrahend,ring)   
   }
   
   fn mont_add(self, addend: Self, ring: Self) -> Self{
      unimplemented!()
   }
   
   fn mont_prod(self, multiplicand: Self, inv: Self, ring: Self) -> Self{
        machine_prime::mont_prod_128(self,multiplicand,ring,inv)
   }
  
   fn mont_sqr(self, inv: Self, ring: Self) -> Self{
        machine_prime::mont_sqr_128(self,ring,inv)
   }
  
   //fn mont_pow(self,one: Self,p: Self, inv: Self, ring: Self) -> Self{
  fn mont_pow(self,mut one: Self,mut p: Self, inv: Self, ring: Self) -> Self{
    machine_prime::mont_pow_128(self,one,p,ring,inv)
   }
   
   fn odd_pow(self, p: Self, n: Self) -> Self{
      
    let one = n.n_identity();
    let base = self.to_mont(n);
    if base < 2{
       return base;
    }
    let npi = n.inv_2();
    base.mont_pow(one,p,npi,n).to_z(npi,n)
   }
  
   fn even_pow(self, mut pow: Self,n: Self) -> Self{
       let mut z = 1;
       let mut base = self;
       while pow > 1 {
           if pow & 1 == 0 {
               base = base.wrapping_mul(base)&n;
                pow >>= 1;
            } else {
                 z = base.wrapping_mul(z)&n;
              base = base.wrapping_mul(base)&n;
               pow = (pow - 1) >> 1
        }
    }
      base.wrapping_mul(z)&n
      }
}

