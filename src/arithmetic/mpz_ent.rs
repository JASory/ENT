use crate::traits::NumberTheory;

use crate::arithmetic::inlineops::*;
use crate::arithmetic::mpz::Mpz;
use crate::arithmetic::sign::Sign;
use crate::arithmetic::sliceops::*;
use crate::data::primes::PRIMELIST;

use crate::data::primes::MERSENNE_LIST;

use std::cmp::Ordering;

impl NumberTheory for Mpz {
    fn rng() -> Self {
        Mpz::unchecked_new(Sign::Positive, vec![rng_64(), rng_64(), rng_64(), rng_64()])
    }

    fn euclidean_div(&self, other: &Self) -> (Self, Self) {
        let (mut quo, mut rem) = self.ref_euclidean(other);
        if self.sign == Sign::Negative && other.sign == Sign::Negative {
            rem.neg();
        } else if self.sign == Sign::Positive && other.sign == Sign::Negative {
            quo.neg();
        } else if self.sign == Sign::Negative && other.sign == Sign::Positive {
            quo.neg();
            rem.neg();
        }

        quo.fix_zero();
        rem.fix_zero();
        (quo, rem)
    }

    fn quadratic_residue(&self, n: &Self) -> Self {
        if n == &Mpz::zero() {
            return Mpz::zero();
        }

        let mut p = self.clone();

        if p.sign == Sign::Negative {
            p = p.add_modinv(n);
        }

        if n.is_power_of_two() {
            let reducer = n.ref_subtraction(&Mpz::one());
            return p.ref_product(&p).and(&reducer);
        }

        p.ref_product(&p).ref_euclidean(n).1
    }

    fn product_residue(&self, other: &Self, n: &Self) -> Self {
        if n == &Mpz::zero() {
            return Mpz::zero();
        }

        let mut p = self.clone();
        let mut q = other.clone();

        if p.sign == Sign::Negative {
            p = p.add_modinv(n);
        }

        if q.sign == Sign::Negative {
            q = q.add_modinv(n)
        }

        if n.is_power_of_two() {
            let reducer = n.ref_subtraction(&Mpz::one());
            return p.ref_product(&q).and(&reducer);
        }

        p.ref_product(&q).ref_euclidean(n).1
    }

    fn exp_residue(&self, y: &Self, n: &Self) -> Self {
        if n == &Mpz::zero() {
            match y.to_u64() {
                Some(x) => return self.pow(x),
                None => panic!("Incomputably large result"),
            }
        }

        let mut p = self.clone();

        if p.sign == Sign::Negative {
            p = p.add_modinv(n);
        }

        p.u_mod_pow(y, n)
    }

    fn checked_exp_residue(&self, y: &Self, n: &Self) -> Option<Self> {
        if n == &Mpz::zero() {
            match y.to_u64() {
                Some(x) => return Some(self.pow(x)),
                None => return None,
            }
        }

        let mut p = self.clone();

        if p.sign == Sign::Negative {
            p = p.add_modinv(n);
        }

        Some(p.u_mod_pow(y, n))
    }

    fn is_sprp(&self, base: &Self) -> bool {
        self.sprp(base.clone())
    }

    #[cfg(not(feature = "parallel"))]
    fn is_prime(&self) -> bool {
        if self.len() < 3 {
            return self.to_u128().unwrap().is_prime();
        }

        if self.is_fermat() {
            return false;
        }

        match self.is_mersenne() {
            Some(x) => {
                if MERSENNE_LIST.contains(&(x as u32)) {
                    return true;
                } else if x < 57885161 {
                    return false;
                } else {
                    // actually incomputable but left for correctness
                    return self.llt(x);
                }
            }
            None => (),
        }

        if !self.trial_div() {
            return false;
        }

        if !self.is_sprp(&Mpz::from_u64(2)) {
            return false;
        }

        if !self.form_check() {
            return false;
        }

        if !self.jacobi_check_mpz() {
            return false;
        }

        if !self.weighted_sprp() {
            return false;
        }
        // if !self.lucas(){
         //return false
      // }
        true
    }
    #[cfg(feature = "parallel")]
    fn is_prime(&self) -> bool {
        if self.len() < 3 {
            return self.to_u128().unwrap().is_prime();
        }
        // if !self.trial_div(){
        //    return false
        // }
        if self.is_fermat() {
            return false;
        }

        match self.is_mersenne() {
            Some(x) => {
                if MERSENNE_LIST.contains(&(x as u32)) {
                    return true;
                } else if x < 57885161 {
                    return false;
                } else {
                    // actually incomputable but left for correctness
                    return self.llt(x);
                }
            }
            None => (),
        }

        if !self.trial_div() {
            return false;
        }
        // println!("Passed the trial div");
        if self.len() < 9 {
            // if less than 2^6400 then serial
            if !self.is_sprp(&Mpz::from_u64(2)) {
                return false;
            }
            if !self.form_check() {
                return false;
            }
            if !self.jacobi_check_mpz() {
                return false;
            }
            if !self.weighted_sprp() {
                return false;
            }
            return true;
        } else {
            // parallelize the checks which roughly cost the same
            // println!("Parallel");
            let a = self.clone();
            let b = self.clone();
            let c = self.clone();
            let two = std::thread::spawn(move || a.is_sprp(&Mpz::from_u64(2)));
            let jacobi = std::thread::spawn(move || b.jacobi_check_mpz());
            let rand = std::thread::spawn(move || c.weighted_sprp());

            if !self.form_check() {
                return false;
            }

            return two.join().unwrap() & jacobi.join().unwrap() & rand.join().unwrap();
        }
    }

    fn prime_proof(&self) -> (bool, Vec<Self>) {
        let x_minus = self.ref_subtraction(&Mpz::one());
        let fctrs = x_minus
            .factor()
            .iter()
            .step_by(2)
            .map(|y| y.clone())
            .collect::<Vec<Self>>();
        let mut certificate = vec![Mpz::from_u64(2)];

        certificate.extend_from_slice(&fctrs[..]);

        loop {
            // loops until it has either been shown to be prime or composite

            let mut witness = Mpz::rand(self.len(), u64::rng)
                .ref_euclidean(&x_minus.ref_subtraction(&Mpz::one()))
                .1
                .ref_addition(&Mpz::from_u64(2)); //Self::rng()%(x-2) + 2;

            'witness: loop {
                if witness.euclid_gcd(self) == Mpz::one() {
                    break 'witness;
                }
                witness.successor();
            }

            if witness.exp_residue(&x_minus, self) != Mpz::one() {
                // If any witness says it's composite then it is
                certificate[0] = witness;

                return (false, certificate);
            }

            'inner: for (idx, j) in fctrs.iter().enumerate() {
                if witness.exp_residue(&x_minus.ref_euclidean(j).0, self) == Mpz::one() {
                    break 'inner;
                }
                if idx == fctrs.len() - 1 {
                    certificate[0] = witness;

                    return (true, certificate);
                }
            }
        }
    }
    //#[cfg(not(feature="parallel"))]
    fn prime_list(&self, sup: &Self) -> Vec<Self> {
        //currently only supports less than, and not inclusive

        let mut delta = self.ref_subtraction(sup);
        delta.successor();
        match delta.to_u64() {
            Some(x) => {
                let mut primevector = vec![];
                let mut p = self.clone();
                for _ in 0..x {
                    if p.is_prime() {
                        primevector.push(p.clone())
                    }
                    p.successor();
                }
                primevector
            }
            None => panic!("Incomputably large interval"),
        }
        //  return Vec::new();
    }
    /*
         #[cfg(feature="parallel")]
        fn prime_list(&self, sup: &Self) -> Vec<Self> {

        }
    */
    fn nth_prime(&self) -> Option<Self> {
        let mut count = Mpz::zero();
        let mut start = Mpz::zero();
        loop {
            start.successor();
            if start.is_prime() {
                count.successor()
            }
            if count.u_cmp(self) == Ordering::Equal {
                return Some(start);
            }
        }
    }

    fn pi(&self) -> Self {
        let mut count = 0u64;
        let mut start = Mpz::zero();
        loop {
            start.successor();
            if start.is_prime() {
                count += 1;
            }
            if start.u_cmp(self) == Ordering::Equal {
                return Mpz::from_u64(count);
            }
        }
    }

    fn prime_gen(x: u32) -> Mpz {
        if x < 128 {
            return Mpz::from_u128(u128::prime_gen(x));
        }
        let mut form = Mpz::one().shl(x as usize - 1);
        let bitlength = form.ref_subtraction(&Mpz::one());
        form.successor();
        let len = (x / 64) as usize;
        loop {
            let mut k = Mpz::rand(len, u64::rng);
            k.mut_and(&bitlength);
            k.mut_or(&form);
            if k.is_prime() {
                return k;
            }
        }
    }

    fn euclid_gcd(&self, other: &Self) -> Self {
        let mut a = self.clone();
        let mut b = other.clone();

        //  let mut t = Mpz::zero();
        while b != Mpz::zero() {
            let t = b.clone();

            b = a.ref_euclidean(&b).1;

            b.normalize();
            a = t;
        }
        a
    }
    // Z variant
    fn extended_gcd(&self, other: &Self) -> (Self, Self, Self) {
        let mut gcd = self.clone();
        let mut new_r = other.clone();
        let mut bezout_1 = Mpz::one();
        let mut new_s = Mpz::zero();
        let mut bezout_2 = Mpz::zero();
        let mut new_t = Mpz::one();

        while new_r != Mpz::zero() {
            let (quo, _rem) = gcd.euclidean_div(&new_r);
            let mut temp = new_r.clone();
            new_r = gcd.ref_subtraction(&quo.ref_product(&temp));
            gcd = temp.clone();

            temp = new_s.clone();
            new_s = bezout_1.ref_subtraction(&quo.ref_product(&temp));
            bezout_1 = temp.clone();

            temp = new_t.clone();
            new_t = bezout_2.ref_subtraction(&quo.ref_product(&temp));
            bezout_2 = temp.clone();
        }
        (gcd, bezout_1, bezout_2)
    }

    fn lcm(&self, other: &Self) -> Self {
        let cf = self.euclid_gcd(other);
        self.ref_product(other).euclidean_div(&cf).1
    }

    fn checked_lcm(&self, other: &Self) -> Option<Self> {
        Some(self.lcm(other))
    }

    //115792089237316195423570985008687907853269984665640564039457584007913129639935
    fn factor(&self) -> Vec<Self> {
        let mut n = self.clone();
        let mut factors: Vec<Self> = vec![];
        let twofactor = n.trailing_zeros();

        if twofactor > 0 {
            n.mut_shr(twofactor as usize);
            factors.push(Mpz::from_u64(2u64));
            factors.push(Mpz::from_u64(twofactor));
        }

        for i in PRIMELIST[1..].iter() {
            // skips two as it has already been eliminated
            let (quo, rem) = n.word_div(*i as u64);

            if rem == 0 {
                let mut count = 1u64;
                factors.push(Mpz::from_u64(*i as u64));
                n = quo;

                'inner: loop {
                    let (inner_quo, inner_rem) = n.word_div(*i as u64);

                    if inner_rem != 0 {
                        break 'inner;
                    }
                    n = inner_quo;
                    count += 1;
                }

                factors.push(Mpz::from_u64(count));
            }
        }
        n.normalize();

        if n == Mpz::one() {
            return factors;
        }

        if n.sprp_check(5) {
            // stops if prime
            factors.push(n.clone());
            factors.push(Mpz::one());
            return factors;
        }

        'outer: while n != Mpz::one() {
            let k = n.rho_mpz();

            factors.push(k.clone());
            let mut count = 0u64;

            'secinner: loop {
                let (inner_quo, inner_rem) = n.ref_euclidean(&k);

                if inner_rem != Mpz::zero() {
                    break 'secinner;
                }
                n = inner_quo;
                n.normalize(); // remove ?
                count += 1;
            }
            factors.push(Mpz::from_u64(count));

            if n.sprp_check(5) {
                // stops if  n is prime
                factors.push(n);
                factors.push(Mpz::one());
                break 'outer;
            }
        }

        factors
    }

    /// Returns the integer component of sqrt(x)
    fn sqrt(&self) -> (Self, Self) {
        if self.len() < 3 {
            let k = Mpz::from_u64(self.to_u128().unwrap().sqrt().0 as u64);
            if !self.is_positive() {
                return (k, Mpz::one());
            }
            return (k, Mpz::zero());
        }
        let mut est = self.shr(((self.bit_length() / 2) - 1) as usize).abs();

        loop {
            let s = est.clone();
            let t = s.ref_addition(&self.euclidean_div(&s).0);
            est = t.shr(1);
            remove_lead_zeros(&mut est.limbs);
            if est.u_cmp(&s) == Ordering::Greater || est.u_cmp(&s) == Ordering::Equal {
                if self.sign == Sign::Negative {
                    return (s, Mpz::one());
                }
                return (s, Mpz::zero());
            }
        }
    }

    fn nth_root(&self, n: &Self) -> (Self, Self) {
        let y = n.to_u64().unwrap();
        let shift = ((self.bit_length() / y) - 1) * (y - 1);
        let mut est = self.shr(shift as usize).abs();
        let scalar = Mpz::from_u64(y - 1);
        let ymp = Mpz::from_u64(y);

        loop {
            let s = est.clone();
            let t = s
                .ref_product(&scalar)
                .ref_addition(&self.euclidean_div(&s.pow(y - 1)).0);
            est = t.euclidean_div(&ymp).0;
            if est.u_cmp(&s) == Ordering::Greater || est.u_cmp(&s) == Ordering::Equal {
                if self.sign == Sign::Negative && self.is_even() {
                    return (s, Mpz::one());
                }
                return (s, Mpz::zero());
            }
        }
    }

    fn radical(&self) -> Self {
        let mut rad = Mpz::one();

        for i in self.factor().iter().step_by(2) {
            rad = rad.ref_product(i)
        }

        rad
    }

    fn k_free(&self, x: &Self) -> bool {
        for i in self.factor()[1..].iter().step_by(2) {
            if i == x {
                return false;
            }
        }
        true
    }

    fn euler_totient(&self) -> Self {
        let mut factors = self.factor();

        let mut denominator = Mpz::one();
        let mut numerator = Mpz::one();

        for i in factors.iter().step_by(2) {
            denominator = denominator.ref_product(i)
        }
        for i in factors.iter_mut().step_by(2) {
            sub_slice(&mut i.limbs[..], &[1]);
            numerator = numerator.ref_product(i)
        }

        (self.ref_euclidean(&denominator).0).ref_product(&numerator)
    }

    fn jordan_totient(&self, k: &Self) -> Option<Self> {
        let fctr = self.factor();

        let mut denom = Mpz::one();
        let mut numer = Mpz::one();
        let negone = Mpz::unchecked_new(Sign::Negative, vec![1]);
        for i in fctr.iter().step_by(2) {
            let pow = i.pow(k.to_u64().unwrap());

            denom = denom.ref_product(&pow);

            numer = numer.ref_product(&pow.ref_addition(&negone));
        }

        Some(
            numer
                .ref_product(&self.pow(k.to_u64().unwrap()))
                .ref_euclidean(&denom)
                .0,
        )
    }

    fn dedekind_psi(&self, k: &Self) -> Option<Self> {
        let k2 = k.shl(1);

        match self.jordan_totient(&k2) {
            Some(y) => Some(y.ref_euclidean(&self.jordan_totient(&k).unwrap()).1),
            None => None,
        }
    }

    fn legendre(&self, p: &Self) -> i8 {
        let mut p_minus = p.clone();
        p_minus.set_sign(Sign::Positive);
        sub_slice(&mut p_minus.limbs[..], &[1]);
        let pow = p_minus.ref_euclidean(&Mpz::from_u64(2)).0;
        let k = self.exp_residue(&pow, p);
        if k == Mpz::one() {
            return 1;
        };
        if k == p_minus {
            return -1;
        };
        0i8
    }

    fn checked_legendre(&self, p: &Self) -> Option<i8> {
        let mut plus = p.clone();
        plus.set_sign(Sign::Positive);
        if plus == Mpz::from_u64(2) {
            return None;
        }
        match plus.is_prime() {
            true => Some(self.legendre(&plus)),
            false => None,
        }
    }

    fn liouville(&self) -> i8 {
        let mut primeomega = Mpz::zero();

        for i in self.factor()[1..].iter().step_by(2) {
            primeomega.mut_addition(i.clone());
        }
        if primeomega.is_even() {
            return 1;
        }
        return -1;
    }

    fn mangoldt(&self) -> f64 {
        let fctr = self.factor();
        if fctr.len() != 2 {
            return 0f64;
        }
        return fctr[0].ln();
    }

    fn jacobi(&self, k: &Self) -> i8 {
        let mut n = self.clone();
        let mut p = k.clone();
        let mut t = 1i8;
        n = n.euclidean_div(&p).1;

        while n != Mpz::zero() {
            let zeros = n.trailing_zeros();
            n.mut_shr(zeros as usize);

            if (p.congruence_u64(8, 3) || p.congruence_u64(8, 5)) && (zeros % 2 == 1) {
                t = -t
            }

            std::mem::swap(&mut n, &mut p);

            if n.congruence_u64(4, 3) && p.congruence_u64(4, 3) {
                t = -t;
            }

            n = n.euclidean_div(&p).1;
        }

        if p == Mpz::one() {
            t
        } else {
            0
        }
    }

    fn checked_jacobi(&self, k: &Self) -> Option<i8> {
        if k.sign == Sign::Positive && k != &Mpz::zero() && !k.is_even() {
            return Some(self.jacobi(k));
        }
        None
    }

    fn smooth(&self) -> Self {
        let k = self.factor();
        k[k.len() - 2].clone()
    }

    fn is_smooth(&self, b: &Self) -> bool {
        let mut k = b.clone();
        k.set_sign(Sign::Positive);
        matches!(self.smooth().u_cmp(&k), Ordering::Less)
    }
}
