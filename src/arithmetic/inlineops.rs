/*

 PRNG

*/
/*

  Precursor functions

*/
#[inline(always)]
pub(crate) const fn split(x: u128) -> (u64, u64) {
    unsafe { std::mem::transmute::<u128, (u64, u64)>(x) }
}
#[inline(always)]
pub(crate) const fn fuse(hi: u64, lo: u64) -> u128 {
    unsafe { std::mem::transmute::<(u64, u64), u128>((lo, hi)) }
}

#[inline]
pub(crate) fn adc(carry: u8, x: u64, y: u64, output: &mut u64) -> u8 {
      #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
      {
       unsafe { core::arch::x86_64::_addcarry_u64(carry, x, y, output) }
      }
     #[cfg(not(any(target_arch = "x86",target_arch="x86_64")))]
     {
       let x128 = x as u128;
       let y128 = y as u128;
       let sum = x128+y128 + carry as u128;
       *output = sum as u64;
       (sum>>64) as u8
     }
}

#[inline]
pub(crate) fn sbb(carry: u8, x: u64, y: u64, output: &mut u64) -> u8 {
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
      {
    unsafe { core::arch::x86_64::_subborrow_u64(carry, x, y, output) }
      }
        #[cfg(not(any(target_arch = "x86",target_arch="x86_64")))]
        {
          let (interim,flag) = x.overflowing_sub(carry.into());
          let (res,flag2) = interim.overflowing_sub(y);
          *output = res;
          if flag || flag2{
             return 1u8 
          }
          0u8
        }
}

#[inline]
pub(crate) fn carry_mul(carry: u64, x: u64, y: u64, output: &mut u64) -> u64 {
    let product = (x as u128 * y as u128) + carry as u128;
    let (jtmp, car) = split(product);
    *output = jtmp;
    car
}

#[inline]
pub(crate) fn mul_acc(carry: u64, x: u64, y: u64, acc: &mut u128) -> u64 {
    *acc += carry as u128;
    *acc += x as u128 * y as u128;
    let lo = *acc as u64;
    *acc >>= 64;
    lo
}

#[inline]
pub(crate) fn carry_div(carry: u64, x: u64, y: u64, output: &mut u64) -> u64 {
    let num = fuse(carry, x);
    let factor = num / y as u128;
    *output = factor as u64;
    (num % y as u128) as u64
}

#[inline]
pub(crate) fn carry_mod(carry: u64, x: u64, y: u64) -> u64 {
    let num = fuse(carry, x);
    (num % y as u128) as u64
}

#[inline]
pub(crate) fn wide_div(upper: u64, lower: u64, divisor: u64) -> (u64, u64) {
    let interim = fuse(upper, lower);
    (
        (interim / divisor as u128) as u64,
        (interim % divisor as u128) as u64,
    )
}

pub(crate) fn divide3by2(ahi: u64, amid: u64, alo: u64, bhi: u64, blo: u64) -> u64 {
    let (mut q0, mut r) = if ahi < bhi {
        let (q0, r) = wide_div(ahi, amid, bhi);
        (q0, r as u128)
    } else {
        (u64::MAX, ahi as u128 + amid as u128)
    };

    while r <= u64::MAX as u128
        && unsafe { std::mem::transmute::<(u64, u64), u128>((alo, r as u64)) }
            < q0 as u128 * blo as u128
    {
        q0 -= 1;
        r += blo as u128;
    }
    q0
}
/*
#[inline]
pub(crate) fn carry_shl(carry: u64, x: u64, places: u32, output: &mut u64) -> u64 {
            *output = (x.overflowing_shl(places).0) | carry;
              unsafe { core::arch::x86_64::_bextr_u64(x, 64 - places, 64) }
}
*/
#[inline]
pub(crate) fn carry_shl(carry: u64, x: u64, places: u32, output: &mut u64) -> u64{
            *output = (x.overflowing_shl(places).0) | carry;
            if places == 0{
              return 0
            }
            x.wrapping_shr(64-places)
}

/*
#[inline]
pub(crate) fn carry_shr(carry: u64, x: u64, places: u32, output: &mut u64) -> u64 {
    *output = (x >> places) | carry;
    if places == 0 {
        return 0;
    }
    unsafe { core::arch::x86_64::_bextr_u64(x, 0, places) << (64 - places) }
}
*/
#[inline]
pub(crate) fn carry_shr(carry: u64, x: u64, places: u32, output: &mut u64) -> u64 {
    *output = (x >> places) | carry;
    if places == 0 {
        return 0;
    }
    x.overflowing_shl(64 - places).0 
}



/*
  Truncation function
*/

pub(crate) fn remove_lead_zeros(x: &mut Vec<u64>) {
    if let Some(&0) = x.last() {
        let len = x.iter().rposition(|&d| d != 0).map_or(0, |i| i + 1);
        x.truncate(len);
    }
    if x.len() < x.capacity() / 4 {
        x.shrink_to_fit();
    }
}
