use crate::ntrait::NumberTheory;
use crate::result::NTResult;
use crate::structs::{Certificate,Factorization};

use crate::arithmetic::inlineops::*;
use crate::arithmetic::mpz::Mpz;
use crate::arithmetic::sign::Sign;
use crate::arithmetic::sliceops::*;
use crate::data::primes::PRIMELIST;

use crate::data::primes::MERSENNE_LIST;

use std::cmp::Ordering;

impl NumberTheory for Mpz {

    fn is_unit(&self) -> bool{
	if self.limbs.len() == 1 && self.limbs[0] == 1{
	  return true;
	}
	false
    }

    fn rng() -> Self {
        Mpz::unchecked_new(Sign::Positive, vec![u64::rng(), u64::rng(), u64::rng(), u64::rng()])
    }
    
    fn residue(&self, ring: Self) -> Self{
      if ring.clone() == Mpz::zero(){
        return self.clone()
      } 
      let rem = self.ref_euclidean(&ring).1;
      if !self.is_positive(){
        return self.ref_subtraction(&rem)
      }
      rem
    }

    fn euclidean_div(&self, other: Self) -> (Self, Self) {
        let (mut quo, mut rem) = self.ref_euclidean(&other);
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

    fn quadratic_residue(&self, n: Self) -> Self {
        if n == Mpz::zero() {
            return Mpz::zero();
        }

        let mut p = self.clone();

        if p.sign == Sign::Negative {
            p = p.add_modinv(&n);
        }

        if n.is_power_of_two() {
            let reducer = n.ref_subtraction(&Mpz::one());
            return p.ref_product(&p).and(&reducer);
        }

        p.ref_product(&p).ref_euclidean(&n).1
    }
    
    fn checked_quadratic_residue(&self, n: Self) -> NTResult<Self>{
        NTResult::Eval(self.quadratic_residue(n))
    }

    fn product_residue(&self, other: Self, n: Self) -> Self {
        if n == Mpz::zero() {
            return Mpz::zero();
        }

        let mut p = self.clone();
        let mut q = other.clone();

        if p.sign == Sign::Negative {
            p = p.add_modinv(&n);
        }

        if q.sign == Sign::Negative {
            q = q.add_modinv(&n)
        }

        if n.is_power_of_two() {
            let reducer = n.ref_subtraction(&Mpz::one());
            return p.ref_product(&q).and(&reducer);
        }

        p.ref_product(&q).ref_euclidean(&n).1
    }
    
    fn checked_product_residue(&self, other: Self, n: Self) -> NTResult<Self>{
       NTResult::Eval(self.product_residue(other,n))
    }
    

    fn exp_residue(&self, y: Self, n: Self) -> Self {
        if n == Mpz::zero() {
            match i64::try_from(y) {
                Ok(x) => {
                  if x < 1i64{
                    panic!("No inverse")
                  }
                return self.pow(x as u64)
                },
                Err(_) => panic!("Incomputably large result"),
            }
        }

        let mut p = self.clone();

        if p.sign == Sign::Negative {
            p = p.add_modinv(&n);
        }

        p.u_mod_pow(&y, &n)
    }

    fn checked_exp_residue(&self, y: Self, n: Self) -> NTResult<Self> {
        if n == Mpz::zero() {
            match i64::try_from(y) {
                Ok(x) => {
                  if x < 0i64{
                    return NTResult::DNE
                  }
                   return NTResult::Eval(self.pow(x as u64)) 
                },
                Err(_) => return NTResult::CompExceeded,
            }
        }

        let mut p = self.clone();

        if p.sign == Sign::Negative {
            p = p.add_modinv(&n);
        }
        
        if !y.is_positive(){
           let (gcd,base,_) = self.eea(n.clone());
           if !gcd.is_unit(){
             return NTResult::DNE
           }
           return NTResult::Eval(base.u_mod_pow(&y,&n))
        }

        NTResult::Eval(p.u_mod_pow(&y, &n))
    }

    fn fermat(&self, base: Self) -> bool{
       base.exp_residue(self.ref_subtraction(&Mpz::one()),self.clone()).is_unit()
    }
    
    
    fn strong_fermat(&self, base: Self) -> bool {
       if self.is_even(){
          return self.fermat(&base);
       }
        self.sprp(base.clone())
    }

    fn is_prime(&self) -> bool {
        if self.len() < 3 {
            return self.to_u128().unwrap().is_prime();
        }
        
        if !self.trial_div() {
            return false;
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
        
    #[cfg(not(feature = "parallel"))]
      {
        if !self.sprp(Mpz::two()) {
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
        true
       } // end single-thread block
           #[cfg(feature = "parallel")]
       {
               if self.len() < 9 {
            // if less than 2^6400 then serial
            if !self.sprp(&Mpz::two()) {
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
            let a = self.clone();
            let b = self.clone();
            let c = self.clone();
            let two = std::thread::spawn(move || a.is_sprp(&Mpz::two()));
            let jacobi = std::thread::spawn(move || b.jacobi_check_mpz());
            let rand = std::thread::spawn(move || c.weighted_sprp());

            if !self.form_check() {
                return false;
            }

            return two.join().unwrap() & jacobi.join().unwrap() & rand.join().unwrap();
        }
       } // end parallel block
    }
    
    fn prime_proof(&self) -> Certificate<Self>{
    
        if self.clone() == Mpz::two() {
            return Certificate::new(self.clone(),Mpz::from(3u64),vec![]);
        }

        let x_minus = self.ref_subtraction(&Mpz::one());
        let fctrs = x_minus.factor().unwrap()
            .factor_iter().map(|x| x_minus.ref_euclidean(x).0)
            .collect::<Vec<Self>>();
         //println!("{:?}",fctrs);
        loop {
            // loops until we get a 

            let mut witness = Mpz::from(u64::rng()).ref_euclidean(&self.ref_subtraction(&Mpz::two())).1.ref_addition(&Mpz::two());
            //println!("witness {}",witness);    
            'witness: loop {

                if witness.coprime(self.clone()) {
                	  break 'witness;
		        }
    
                witness.successor();
            }
              
	   if !self.strong_fermat(witness.clone()){
		return Certificate::new(self.clone(),witness.clone(),fctrs)	
	   }

           
            'inner: for (idx, i) in fctrs.iter().enumerate() {
                if witness.exp_residue(i.clone(), self.clone()).is_unit(){
                    break 'inner;
                }
                if idx == fctrs.len() - 1 {
                    return Certificate::new(self.clone(),witness.clone(),fctrs);
                }
            }
	  }
    
    
    }
    
    //#[cfg(not(feature="parallel"))]
    fn prime_list(&self, sup: Self) -> Vec<Self> {
        //currently only supports less than, and not inclusive

        let mut delta = self.ref_subtraction(&sup);
        delta.successor();
        match u64::try_from(delta) {
            Ok(x) => {
                let mut primevector = vec![];
                let mut p = self.clone();
              //  let primevec = 17000u64.prime_list(&30_000);
                for _ in 0..x{
                    if p.is_prime() {
                        primevector.push(p.clone())
                    }
                    p.successor();
                }
                primevector
            }
            Err(_) => panic!("Incomputably large interval"),
        }
        //  return Vec::new();
    }
    
    fn nth_prime(&self) -> NTResult<Self> {
       
        let mut count = Mpz::zero();
        let mut start = Mpz::zero();
        if *self == Mpz::zero(){
          return NTResult::DNE
        }
        loop {
            start.successor();
            if start.is_prime() {
                count.successor()
            }
            if count.u_cmp(self) == Ordering::Equal {
                return NTResult::Eval(start);
            }
        }
    }


    fn pi(&self) -> Self {
      match self.to_u128(){
         Some(x) => return Mpz::from(x.pi()),
         None => (),
      }
        let mut count = 0u64;
        let mut start = Mpz::zero();
        loop {
            start.successor();
            if start.is_prime() {
                count += 1;
            }
            if start.u_cmp(self) == Ordering::Equal {
                return Mpz::from(count);
            }
        }
    }


    fn prime_gen(x: u32) -> NTResult<Mpz> {
        if x < 128 {
         return u128::prime_gen(x).map(Mpz::from)  
        }
        
        let mut form = Mpz::one().shl(x as usize-1);

        let bitlength = form.ref_subtraction(&Mpz::one());

        form.successor();

        loop {
            let mut k = Mpz::rand(x as usize+1);
        
            k.mut_and(&bitlength);

            k.mut_or(&form);
            assert_eq!(k.bit_length() as u32,x);
            if k.is_prime() {
                return NTResult::Eval(k);
            }
        }
    } 

    fn gcd(&self, other: Self) -> Self {
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
    
    // Z variant
    fn extended_gcd(&self, other: Self) -> (Self, Self, Self) {
        let mut gcd = self.clone();
        let mut new_r = other.clone();
        let mut bezout_1 = Mpz::one();
        let mut new_s = Mpz::zero();
        let mut bezout_2 = Mpz::zero();
        let mut new_t = Mpz::one();

        while new_r != Mpz::zero() {
            let (quo, _rem) = gcd.euclidean_div(new_r.clone());
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

    fn lcm(&self, other: Self) -> NTResult<Self> {
        if self == &Mpz::zero() && other  == Mpz::zero(){
          return NTResult::Eval(Mpz::zero())
        }
        
        let cf = self.gcd(other.clone());
        NTResult::Eval(self.euclidean_div(cf).0.ref_product(&other))
    }

// FIXME 
    fn factor(&self) -> NTResult<Factorization<Self>> {
        let mut n = self.clone();
        let mut factors: Factorization<Self> = Factorization::new();
        let twofactor = n.trailing_zeros();

        if twofactor > 0 {
            n.mut_shr(twofactor as usize);
            factors.add_factor(Mpz::from(2u64));
            factors.add_power(twofactor);
        }

        for i in PRIMELIST[1..].iter() {
            // skips two as it has already been eliminated
            let (quo, rem) = n.word_div(*i as u64);

            if rem == 0 {
                let mut count = 1u64;
                factors.add_factor(Mpz::from(*i as u64));
                n = quo;

                'inner: loop {
                    let (inner_quo, inner_rem) = n.word_div(*i as u64);

                    if inner_rem != 0 {
                        break 'inner;
                    }
                    n = inner_quo;
                    count += 1;
                }

                factors.add_power(count);
            }
        }
        n.normalize();

        if n.is_unit(){
             //       println!("Final {}",factors);
            return NTResult::Eval(factors);
        }

        if n.sprp_check(5) {
            // stops if prime
            factors.add_factor(n.clone());
            factors.add_power(1);
            return NTResult::Eval(factors);
        }

        'outer: while !n.is_unit() {
            let k = n.rho_mpz();

            factors.add_factor(k.clone());
            let mut count = 0u64;

            'secinner: loop {
                let (inner_quo, inner_rem) = n.ref_euclidean(&k);

                if !inner_rem.is_zero() {
                    break 'secinner;
                }
                n = inner_quo;
                n.normalize(); // remove ?
                count += 1;
            }
            factors.add_power(count);

            if n.sprp_check(5) {
                // stops if  n is prime
                factors.add_factor(n);
                factors.add_power(1);
                break 'outer;
            }
        }

        NTResult::Eval(factors)
    }
    

    /// Returns the integer component of sqrt(x)
    fn sqrt(&self) -> (Self, Self) {
        if self.len() < 3 {
            let k = Mpz::from(self.to_u128().unwrap().sqrt().0 as u64);
            if !self.is_positive() {
                return (k, Mpz::one());
            }
            return (k, Mpz::zero());
        }
        let mut est = self.shr(((self.bit_length() / 2) - 1) as usize).abs();

        loop {
            let s = est.clone();
            let t = s.ref_addition(&self.euclidean_div(s.clone()).0);
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

    fn nth_root(&self, n: Self) -> (Self, Self) {
        let y = u64::try_from(n).unwrap();
        let shift = ((self.bit_length() / y) - 1) * (y - 1);
        let mut est = self.shr(shift as usize).abs();
        let scalar = Mpz::from(y - 1);
        let ymp = Mpz::from(y);

        loop {
            let s = est.clone();
            let t = s
                .ref_product(&scalar)
                .ref_addition(&self.euclidean_div(s.pow(y - 1)).0);
            est = t.euclidean_div(ymp.clone()).0;
            if est.u_cmp(&s) == Ordering::Greater || est.u_cmp(&s) == Ordering::Equal {
                if self.sign == Sign::Negative && self.is_even() {
                    return (s, Mpz::one());
                }
                return (s, Mpz::zero());
            }
        }
    }
    // FIXME : Use prime exponents and fix for x < 0
   fn max_exp(&self) -> (Self,Self){
     let mut expo = Mpz::from(self.bit_length());
     let mut flag = true;
     if  self.is_positive(){
       flag = false
     }
      loop{
        let mut base = self.abs().nth_root(expo.clone()).0;
         if base.pow(u64::try_from(expo.clone()).unwrap()) == self.abs(){
         if flag{
           base.neg();
         }
           return(base,expo.clone())
         }
         expo.predecessor();
         if u64::try_from(expo.clone()).unwrap() == 1{
         let mut val = self.clone();
          if flag {
            val.neg()
          }
           return (val,Mpz::one())
         }
      }
    }


    
    fn radical(&self) -> NTResult<Self>{
       match u64::try_from(self.clone()){
        Ok(x) => x.radical().map(Mpz::from),
        Err(_) => {
          let mut rad = Mpz::one();
           for i in self.factor().unwrap().factor_iter(){
              rad = rad.ref_product(i)
           }
          NTResult::Eval(rad)
         },
       }
    }

    fn k_free(&self, x: Self) -> bool {
       let x = u64::try_from(x).unwrap();
        for i in self.factor().unwrap().power_iter(){
            if *i == x {
                return false;
            }
        }
        true
    }

    fn euler_totient(&self) -> Self {
    
       match self.to_u128(){
      
         Some(x)=> return Mpz::from(x.euler_totient()),
         None => (),
       }
       
        let mut factors = self.factor().unwrap();

        let mut denominator = Mpz::one();
        let mut numerator = Mpz::one();

        for i in factors.factor_iter() {
            denominator = denominator.ref_product(i)
        }
        for i in factors.factor_iter() {
            let mut fctr = i.clone();
            fctr.predecessor();
            numerator = numerator.ref_product(&fctr);
        }

        (self.ref_euclidean(&denominator).0).ref_product(&numerator)
    }

    fn jordan_totient(&self, k: Self) -> NTResult<Self> {
        if self.clone() == Mpz::zero() || self.is_unit(){
          return NTResult::Eval(self.clone())
        }
        let fctr = self.factor().unwrap();
        let mut denom = Mpz::one();
        let mut numer = Mpz::one();
        let negone = Mpz::unchecked_new(Sign::Negative, vec![1]);
        let value = u64::try_from(k).unwrap();
        for i in fctr.factor_iter() {
            let pow = i.pow(value);

            denom = denom.ref_product(&pow);

            numer = numer.ref_product(&pow.ref_addition(&negone));
        }

        NTResult::Eval(
            numer
                .ref_product(&self.pow(value))
                .ref_euclidean(&denom)
                .0,
        )
    }
    
    fn exponent(&self) -> NTResult<Self>{
       
       match self.to_u128(){
         Some(x)=> return x.exponent().map(Mpz::from),
         None => (),
       }
       
       let fctr = self.factor().unwrap();
       //let base = fctr.iter()..map(|z| z.clone()).collect::<Vec<Self>>();
       let power = fctr.power_iter().cloned().collect::<Vec<u64>>();
       let two = Mpz::two();
       let mut result = Mpz::one();
      for (idx,el) in fctr.factor_iter().enumerate(){
        if el == &two && power[0] > 2{
         let phi = ((el.pow(power[idx]).ref_euclidean(el).0).ref_product(&el.ref_subtraction(&Mpz::one()))).shr(1);
          result = result.lcm(phi).unwrap();
        }
       else{
         let phi = (el.pow(power[idx]).ref_euclidean(el).0).ref_product(&el.ref_subtraction(&Mpz::one()));
         result = result.lcm(phi).unwrap();
       } 
      }
     NTResult::Eval(result)
    }

    fn dedekind_psi(&self, k: Self) -> NTResult<Self> {
     if *self == Mpz::zero(){
         return NTResult::Infinite
       }
        let k2 = k.shl(1);
     self.jordan_totient(k2).map(|y| y.ref_euclidean(&self.jordan_totient(k).unwrap()).0)
    }

    fn legendre(&self, p: Self) -> i8 {
        let mut p_minus = p.clone();
        p_minus.set_sign(Sign::Positive);
        p_minus.predecessor();

        let pow = p_minus.ref_euclidean(&Mpz::two()).0;
        let k = self.exp_residue(pow, p);

        if k == Mpz::one() {
            return 1;
        };
        if k == p_minus {
            return -1;
        };
        0i8
    }

    fn checked_legendre(&self, p: Self) -> NTResult<i8> {
        let mut plus = p.clone();
        plus.set_sign(Sign::Positive);
        if plus == Mpz::two() {
            return NTResult::Undefined;
        }
        match plus.is_prime() {
            true => NTResult::Eval(self.legendre(plus)),
            false => NTResult::Undefined,
        }
    }

    fn liouville(&self) -> i8 {
        match self.to_u128(){
          Some(an) => return an.liouville(),
          None => (),
        }
        
        let mut primeomega = 0;

        for i in self.factor().unwrap().power_iter() {
            primeomega+=i;
        }
        if primeomega&1==0 {
            return 1;
        }
        -1
    }
    
    fn derivative(&self) -> NTResult<Self> {
       if *self == Mpz::zero(){
         return NTResult::Eval(Mpz::zero())
       }
       if self.is_unit(){
        return NTResult::Eval(Mpz::zero())
       }
       let fctr = self.factor().unwrap();
       let mut sum = Mpz::zero();
       
     for i in 0..fctr.base.len() / 2 {
        sum.mut_addition(Mpz::from(fctr.power[i]).ref_product(&self.ref_euclidean(&fctr.base[i]).0))
      }
    NTResult::Eval(sum)
    }

    fn mangoldt(&self) -> f64 {
     match self.to_u128(){
          Some(eburum) => return eburum.mangoldt(),
          None => (),
        }
        let fctr = self.factor().unwrap();
        if fctr.base.len() != 2 {
            return 0f64;
        }
        fctr.base[0].ln()
    }
    
    
    fn mobius(&self) -> i8 {
     match self.to_u128(){
          Some(eburum) => return eburum.mobius(),
          None => (),
        }
      let fctr = self.factor().unwrap();
      /*
      if fctr.len() == 1{ // if only one factor then return -1
         return -1
      }
      */
      for i in fctr.power_iter(){
        if *i  > 1{
         return 0
        }
      }
      
      let fctrsum = fctr.power_iter().sum::<u64>();
      if fctrsum&1 == 1{// if odd number of factors and square free
        return -1
      }
      1
    }

    fn jacobi(&self, k: Self) -> i8 {
        let mut n = self.clone();
        let mut p = k.clone();
        let mut t = 1i8;
        n = n.euclidean_div(k).1;

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

            n = n.euclidean_div(p.clone()).1;
        }

        if p == Mpz::one() {
            t
        } else {
            0
        }
    }

    fn checked_jacobi(&self, k: Self) -> NTResult<i8> {
        if k.sign == Sign::Positive && k != Mpz::zero() && !k.is_even() {
            return NTResult::Eval(self.jacobi(k));
        }
        NTResult::Undefined
    }
    
     fn kronecker(&self, k: Self) -> i8{
     let mut multiplier = 1;
    
    if !k.is_positive() && !self.is_positive(){
      multiplier = -1;
    }
     let x = self.clone();
     if k.is_zero(){
      if x.is_unit(){
         return 1
      }
     return 0
    }
   if k.is_unit(){
      return 1
   }
   let fctr = k.factor().unwrap();
   let mut start = 0;
   let mut res = 1;
   let two = Mpz::two();
   if fctr.base[0] == two{
     start = 1;
     if x.is_even(){
     res = 0;
     }
     else if x.congruence_u64(8,1) || x.congruence_u64(8,7){
      res=1
     }
     else{
       res = (-1i8).pow(fctr.power[1] as u32)
     }
   }
   if fctr.base[0] == two && fctr.base.len() == 2{
     return res*multiplier
   }
   for i in start..fctr.base.len(){
     res*=self.legendre(fctr.base[i].clone()).pow(fctr.power[i] as u32);
   }
   res*multiplier
}
    
    fn smooth(&self) -> NTResult<Self> {
       if *self == Mpz::zero(){
         return NTResult::Infinite
       }
       if self.is_unit(){
        return NTResult::DNE
       }
      NTResult::Eval(self.factor().unwrap().max())
    }
    
    

    fn is_smooth(&self, b: Self) -> bool {
     match self.smooth(){
      NTResult::Infinite => false,
      NTResult::Eval(x) => {
      let mut k = b.clone();
        k.set_sign(Sign::Positive);
        matches!(x.u_cmp(&k), Ordering::Less)
      }, 
      _=> false,
     }
   }
   
   fn ord(&self, n: Self) -> NTResult<Self>{

        let ord_2 = |a: Mpz, p: u64| -> Mpz{
        let modulo = Mpz::one().shl(p as usize);
        let mut b = a.ref_euclidean(&modulo).1;
   
         if b == Mpz::one(){
            return Mpz::one();
         }
      for i in 1..p{
         b = b.quadratic_residue(modulo.clone());
         if b == Mpz::one(){
           return Mpz::one().shl(i as usize);
         }
    }
    return p.into();
    };

// Given ord(a,p)  calculate ord(a,p^n)
    let  pp_ord = |a: Mpz, b: Mpz, p: Mpz, e: u64| -> Mpz{
        let fullpow = p.pow(e);

     	for i in 0..e+1{
          let partial = b.ref_product(&p.pow(i));
          if a.exp_residue(partial.clone(),fullpow.clone()) == Mpz::one(){
             return partial;
           }
        }
    return fullpow;
    };

    let p_ord = |a: Mpz, p: Mpz| -> Mpz{

 	  let mut pminus = p.ref_subtraction(&Mpz::one());  
  	  let fctr = pminus.factor().unwrap();
   	  let mut m = pminus.clone();

 	for i in fctr.pair_iter(){

 	    for _ in 0..*i.1{
              let reduced = m.ref_euclidean(&i.0).0;

	        if a.exp_residue(reduced.clone(),p.clone()) == Mpz::one(){
                   m = reduced;
                }
          	else{
                   break;
          	}
            }
        }
       m
     };


    if !self.coprime(n.clone()){
       return NTResult::DNE;
    }

    let fctr = n.clone().factor().unwrap();
    let mut fullord = Mpz::one();

    for i in fctr.pair_iter(){
     let mut ord : Self;
      if *i.0 == Mpz::two(){
         ord = ord_2(self.clone(),*i.1);
      }
      else{
        ord = p_ord(self.clone(),i.0.clone());
        if *i.1 > 1{
           ord=pp_ord(self.clone(),ord,i.0.clone(),*i.1);
        }
      }
       fullord = fullord.lcm(ord).unwrap(); 
    }
    NTResult::Eval(fullord)
       
   }
   
   fn mul_inverse(&self, n: Self) -> NTResult<Self>{
       unimplemented!()
   }
   
}
