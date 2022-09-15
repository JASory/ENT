/*
   Implementation of Montgomery arithmetic
*/

use crate::data::primes::INV_8;
use crate::primitive::sixteenbytes::*;

// use crate::traits::NumberTheory;

fn mod_inv32(n: u32) -> u32 {
    // inverse of odd n in  2^32
    let mut est = INV_8[((n >> 1) & 0x7F) as usize] as u32;
    est = 2u32.wrapping_sub(est.wrapping_mul(n)).wrapping_mul(est);
    est = 2u32.wrapping_sub(est.wrapping_mul(n)).wrapping_mul(est);
    est = 2u32.wrapping_sub(est.wrapping_mul(n)).wrapping_mul(est);
    est.wrapping_neg()
}

pub(crate) fn mod_inv64(n: u64) -> u64 {
    // inverse of odd n in  2^64
    let mut est = INV_8[((n >> 1) & 0x7F) as usize] as u64;
    est = 2u64.wrapping_sub(est.wrapping_mul(n)).wrapping_mul(est);
    est = 2u64.wrapping_sub(est.wrapping_mul(n)).wrapping_mul(est);
    est = 2u64.wrapping_sub(est.wrapping_mul(n)).wrapping_mul(est);
    est.wrapping_neg()
}

fn mod_inv128(n: u128) -> u128 {
    // inverse of odd n in  2^128
    let mut est = INV_8[((n >> 1) & 0x7F) as usize] as u128;
    est = 2u128.wrapping_sub(est.wrapping_mul(n)).wrapping_mul(est);
    est = 2u128.wrapping_sub(est.wrapping_mul(n)).wrapping_mul(est);
    est = 2u128.wrapping_sub(est.wrapping_mul(n)).wrapping_mul(est);
    est = 2u128.wrapping_sub(est.wrapping_mul(n)).wrapping_mul(est);
    est = 2u128.wrapping_sub(est.wrapping_mul(n)).wrapping_mul(est);
    est.wrapping_neg()
}

fn x32_modn(x: u32, n: u32) -> u32 {
    (((x as u64) << 32) % (n as u64)) as u32
}

pub fn x64_modn(x: u64, n: u64) -> u64 {
    (((x as u128) << 64) % (n as u128)) as u64
}

pub fn x128_modn(x: u128, n: u128) -> u128 {
    u256mod128((x, 0), n)
    //n.wrapping_neg()
    //(u128::MAX%n)+1
}

fn montprod_32(x: u32, y: u32, n: u32, npi: u32) -> u32 {
    let input = x as u64 * y as u64;
    let tm = (input as u32).wrapping_mul(npi);
    let (t, overflow) = input.overflowing_add((tm as u64) * (n as u64));
    let t = (t >> 32) as u32;

    if overflow {
        t + n.wrapping_neg()
    } else if t >= n {
        t - n
    } else {
        t
    }
}

fn montprod_64(x: u64, y: u64, n: u64, npi: u64) -> u64 {
    let input = x as u128 * y as u128;
    let tm = (input as u64).wrapping_mul(npi);
    let (t, overflow) = input.overflowing_add((tm as u128) * (n as u128));
    let t = (t >> 64) as u64;

    if overflow {
        t + n.wrapping_neg()
    } else if t >= n {
        t - n
    } else {
        t
    }
}

fn montprod_128(x: u128, y: u128, n: u128, npi: u128) -> u128 {
    let (phi, plo) = u256prod(x, y);

    let tm = plo.wrapping_mul(npi);

    let (t, overflow) = overflowing_add(u256prod(n, tm), (phi, plo));

    let t = t.0;

    if overflow {
        t + n.wrapping_neg()
    } else if t >= n {
        t - n
    } else {
        t
    }
}

fn mpow_32(x: u32, p: u32, n: u32, npi: u32) -> u32 {
    let mut z = x32_modn(1u32, n);
    let mut base = x32_modn(x, n);
    let mut pow = p;

    while pow > 1 {
        if pow & 1 == 0 {
            base = montprod_32(base, base, n, npi);
            pow >>= 1;
        } else {
            z = montprod_32(base, z, n, npi);
            base = montprod_32(base, base, n, npi);
            pow = (pow - 1) >> 1
        }
    }
    montprod_32(base, z, n, npi)
}

fn mpow_64(x: u64, p: u64, n: u64, npi: u64) -> u64 {
    let mut z = x64_modn(1u64, n);
    let mut base = x64_modn(x, n);
    let mut pow = p;

    while pow > 1 {
        if pow & 1 == 0 {
            base = montprod_64(base, base, n, npi);
            pow >>= 1;
        } else {
            z = montprod_64(base, z, n, npi);
            base = montprod_64(base, base, n, npi);
            pow = (pow - 1) >> 1
        }
    }
    montprod_64(base, z, n, npi)
}

fn mpow_128(x: u128, p: u128, n: u128, npi: u128) -> u128 {
    let mut z = (u128::MAX % n) + 1; //x128_modn(1u128, n);
    let mut base = x128_modn(x, n);
    let mut pow = p;

    while pow > 1 {
        if pow & 1 == 0 {
            base = montprod_128(base, base, n, npi);
            pow >>= 1;
        } else {
            z = montprod_128(base, z, n, npi);
            base = montprod_128(base, base, n, npi);
            pow = (pow - 1) >> 1
        }
    }
    montprod_128(base, z, n, npi)
}

pub(crate) fn sprp_128(p: u128, base: u128) -> bool {
    let p_minus = p - 1;
    let zeroes = p_minus.trailing_zeros();
    let d = p_minus >> zeroes;

    let npi = mod_inv128(p);
    let mut x = mpow_128(base, d, p, npi);
    let one = (u128::MAX % p) + 1; //x128_modn(1, p);
    let oneinv = x128_modn(p_minus, p);

    if x == one || x == oneinv {
        return true;
    }
    for _ in 1..zeroes {
        x = montprod_128(x, x, p, npi);

        if x == oneinv {
            return true;
        }
    }
    false
}

pub(crate) fn sprp_32(p: u32, base: u32) -> bool {
    let p_minus = p - 1;
    let zeroes = p_minus.trailing_zeros();
    let d = p_minus >> zeroes;

    let npi = mod_inv32(p);
    let mut x = mpow_32(base, d, p, npi);
    let one = x32_modn(1, p);
    let oneinv = x32_modn(p_minus, p);

    if x == one || x == oneinv {
        return true;
    }
    for _ in 1..zeroes {
        x = montprod_32(x, x, p, npi);

        if x == oneinv {
            return true;
        }
    }
    false
}

pub(crate) fn sprp_64(p: u64, base: u64) -> bool {
    let p_minus = p - 1;
    let zeroes = p_minus.trailing_zeros();
    let d = p_minus >> zeroes;

    let npi = mod_inv64(p);
    let mut x = mpow_64(base, d, p, npi);
    let one = x64_modn(1, p);
    let oneinv = x64_modn(p_minus, p);
    if x == one || x == oneinv {
        return true;
    }
    for _ in 1..zeroes {
        x = montprod_64(x, x, p, npi);

        if x == oneinv {
            return true;
        }
    }
    false
}

// Fuses two strong fermat tests to minimize total operations. Reference implementation, not used
#[inline]
pub(crate) fn _sprp_64_double(p: u64, base1: u64, base2: u64) -> bool {
    let (mut ef_1, mut ef_2) = (false, false);
    let p_minus = p - 1;
    let twofactor = p_minus.trailing_zeros();
    let mut d = p_minus >> twofactor;

    let npi = mod_inv64(p);
    let one = (u64::MAX % p) + 1;
    let mut z1 = one;
    let mut z2 = one;

    let mut result1 = (((base1 as u128) << 64) % (p as u128)) as u64;
    let mut result2 = (((base2 as u128) << 64) % (p as u128)) as u64;

    let oneinv = (((p_minus as u128) << 64) % (p as u128)) as u64;

    while d > 1 {
        if d & 1 == 0 {
            result1 = montprod_64(result1, result1, p, npi);
            result2 = montprod_64(result2, result2, p, npi);
            d >>= 1;
        } else {
            z1 = montprod_64(z1, result1, p, npi);
            result1 = montprod_64(result1, result1, p, npi);

            z2 = montprod_64(z2, result2, p, npi);
            result2 = montprod_64(result2, result2, p, npi);

            d = (d - 1) >> 1;
        }
    }

    result1 = montprod_64(z1, result1, p, npi);
    result2 = montprod_64(z2, result2, p, npi);

    if result1 == one || result1 == oneinv {
        ef_1 = true;
    }

    if result2 == one || result2 == oneinv {
        ef_1 = true;
    }

    if ef_1 & ef_2 {
        return true;
    }

    for _ in 1..twofactor {
        result1 = montprod_64(result1, result1, p, npi);
        result2 = montprod_64(result2, result2, p, npi);

        if result1 == oneinv {
            ef_1 = true;
        }
        if result2 == oneinv {
            ef_2 = true;
        }
        if ef_1 & ef_2 {
            return true;
        }
    }
    false
}

fn odd_pow32(x: u32, p: u32, n: u32) -> u32 {
    let npi = mod_inv32(n);
    let interim = mpow_32(x, p, n, npi);
    montprod_32(1u32, interim, n, npi)
}

fn odd_pow64(x: u64, p: u64, n: u64) -> u64 {
    let npi = mod_inv64(n);
    let interim = mpow_64(x, p, n, npi);
    montprod_64(1u64, interim, n, npi)
}

fn odd_pow128(x: u128, p: u128, n: u128) -> u128 {
    let npi = mod_inv128(n);
    let interim = mpow_128(x, p, n, npi);
    montprod_128(1u128, interim, n, npi)
}

pub(crate) fn pow_32(x: u32, y: u32, n: u32) -> u32 {
    if n & 1 == 0 {
        let k = n.trailing_zeros() as u64;
        let s = n >> k;

        let reducer = (1 << k) - 1; // A shorthand for arithmetic over Z[2k]

        let k_rem = x.wrapping_pow(y) & reducer;

        let s_rem = odd_pow32(x, y, s);

        let mut s_inv = s;

        for _ in 0..10 {
            // Multiplicative inverse over Z[2k]
            s_inv = 2u32.wrapping_sub(s_inv.wrapping_mul(s)).wrapping_mul(s_inv) & reducer;
        }

        let y = k_rem.wrapping_sub(s_rem).wrapping_mul(s_inv) & reducer;

        s_rem + s * y
    } else {
        odd_pow32(x, y, n)
    }
}

fn even_pow_64(x: u64, y: u64, reducer: u64) -> u64 {
    let mut z = 1u64;
    let mut base = x;

    let mut pow = y;
    if pow == 0 {
        return z;
    }

    while pow > 1 {
        if pow % 2 == 0 {
            base = base.wrapping_mul(base);
            pow >>= 1;
        } else {
            z = base.wrapping_mul(z);
            base = base.wrapping_mul(base);
            pow = (pow - 1) >> 1;
        }
    }

    base.wrapping_mul(z) & reducer
}

fn even_pow_128(x: u128, y: u128, reducer: u128) -> u128 {
    let _result = x;
    let mut z = 1u128;
    let mut base = x;
    // let n = modulus as u128;
    let mut pow = y;
    if pow == 0 {
        return z;
    }

    while pow > 1 {
        if pow % 2 == 0 {
            base = base.wrapping_mul(base);
            pow >>= 1;
        } else {
            z = base.wrapping_mul(z);
            base = base.wrapping_mul(base);
            pow = (pow - 1) >> 1;
        }
    }

    base.wrapping_mul(z) & reducer
}

pub(crate) fn pow_64(x: u64, y: u64, n: u64) -> u64 {
    if n & 1 == 0 {
        let k = n.trailing_zeros() as u64;
        let s = n >> k;

        let reducer = (1 << k) - 1; // A shorthand for arithmetic over Z[2k]

        let k_rem = even_pow_64(x, y, reducer); //x.wrapping_pow(y as u32)&reducer;

        let s_rem = odd_pow64(x, y, s);

        let mut s_inv = s;

        for _ in 0..10 {
            // Multiplicative inverse over Z[2k]
            s_inv = 2u64.wrapping_sub(s_inv.wrapping_mul(s)).wrapping_mul(s_inv) & reducer;
        }

        let y = k_rem.wrapping_sub(s_rem).wrapping_mul(s_inv) & reducer;

        s_rem + s * y
    } else {
        odd_pow64(x, y, n)
    }
}

pub(crate) fn pow_128(x: u128, y: u128, n: u128) -> u128 {
    if n & 1 == 0 {
        let k = n.trailing_zeros() as u128;
        let s = n >> k;

        let reducer = (1 << k) - 1; // A shorthand for arithmetic over Z[2k]

        let k_rem = even_pow_128(x, y, reducer);

        let s_rem = odd_pow128(x, y, s);

        let mut s_inv = s;

        for _ in 0..10 {
            // Multiplicative inverse over Z[2k]
            s_inv = 2u128
                .wrapping_sub(s_inv.wrapping_mul(s))
                .wrapping_mul(s_inv)
                & reducer;
        }

        let y = k_rem.wrapping_sub(s_rem).wrapping_mul(s_inv) & reducer;

        s_rem + s * y
    } else {
        odd_pow128(x, y, n)
    }
}
