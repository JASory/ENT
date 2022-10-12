/*

  Primality

*/

use crate::data::primes::PRIMELIST;

use crate::traits::NumberTheory;
use crate::Mpz;

use crate::arithmetic::sliceops::mod_slice;
use crate::arithmetic::sliceops::sub_slice;


impl Mpz {
    /// Strong Fermat test to a selected base
    pub(crate) fn sprp(&self, base: Self) -> bool {
        let mut p_minus = self.clone();
        let one = Mpz::one();

        sub_slice(&mut p_minus.limbs[..], &one.limbs[..]); //subtract one this will always succeed

        let zeroes = p_minus.trailing_zeros() as usize;

        let d = p_minus.shr(zeroes);
        let mut x = base.u_mod_pow(&d, self);

        if x == Mpz::one() || x == p_minus {
            return true;
        }

        for _ in 0..zeroes - 1 {
            x = x.u_quadratic_residue(self);

            if x == p_minus {
                return true;
            }
        }
        false
    }

    #[cfg(not(feature = "parallel"))]
    /** Performs n random base checks, can be combined with is_prime to strengthen it. As is_prime is equivalent to one random sprp check in the worst case, the function below is 2^-64 in the worst case

     use number_theory::NumberTheory;
     use number_theory::Mpz;

     fn strong_check(x: &Mpz) -> bool{
       if !x.is_prime(){
         return false
       }
       x.sprp_check(31)
     }

    */
    pub fn sprp_check(&self, n: usize) -> bool {
        if self.len() < 2 {
            return self.to_u64().unwrap().is_prime();
        } // if fits into u64, reduce to 64-bit check

        if self.is_fermat() {
            return false;
        }

        for _ in 0..n {
            let z = Mpz::rand((self.bit_length()-1) as usize);

            if !self.sprp(z) {
                return false;
            } //
        }

        true
    }

    #[cfg(feature = "parallel")] // Fix this,likely does not check the correct amount of tests
    /// Checks n bases 
    pub fn sprp_check(&self, k: usize) -> bool {
        if self.len() < 2 {
            return self.to_u64().unwrap().is_prime();
        } // if fits into u64, reduce to 64-bit check

        if self.is_fermat() {
            return false;
        }

       
        let single = |x: Mpz, n: usize| {
            for _ in 0..n {
                let z = Mpz::rand(x.bit_length() as usize-1);
                if !x.sprp(z) {
                    return false;
                }
            }
            return true;
        };

        let threadcount = usize::from(std::thread::available_parallelism().unwrap()) - 2;
        let q = self.clone();
        /*  if k < threadcount{
           let mut threadvec = vec![];
           let mut xclone = vec![];
           for i in 0..k{
             xclone.push(self.clone())
           }
           for i in xclone{
             threadvec.push(std::thread::spawn(move|| {single(i, 1).clone()}));
           }
           for j in threadvec{
             if !j.join().unwrap(){
               return false
             }
           }
           return true
         }
        // else{*/
        let mut threadvec = vec![];
        let mut xclone = vec![];
        for _ in 0..threadcount {
            xclone.push(self.clone())
        }

        for i in xclone {
            threadvec.push(std::thread::spawn(move || single(i, k / threadcount)));
        }
        let tally = single(q, k / threadcount); //threadvec.push(std::thread::spawn(move || {single(q.clone(), k/threadcount)}));

        for j in threadvec {
            if !j.join().unwrap() {
                return false;
            }
        }
        return tally;
    }

    /** Detects if self is a number of various forms and returns the prime evaluation if it is
     One might be tempted to ask why proth numbers are not evaluated, and it's simply the fact that Proth's theorem is actually an extremely
     inefficient  primality test, in the worst case strong fermat tests are twice as efficient as a Proth test and far stronger in practice.
     */
    pub(crate) fn form_check(&self) -> bool {
        // Detects if self is of the form r(k^2 + k) + k + 1

        let sqrt = self.sqrt().0;
        if sqrt.sqr() == self.clone() {
            // detect perfect square
            return false;
        }

        !detect_pseudo_mpz(self)
    }

    pub(crate) fn trial_div(&self) -> bool {
        //weighted trial division

        if self.is_even() {
            return false;
        }

        let rem = mod_slice(&self.limbs[..], 16294579238595022365);

        for i in PRIMELIST[1..16usize].iter() {
            if rem % *i as u64 == 0 {
                return false;
            }
        }

        let mut supremum: usize = 2048;

        if self.len() < 40 {
            supremum = self.len() * 50
        }
        for i in PRIMELIST[17..supremum].iter() {
            if self.congruence_u64(*i as u64, 0) {
                return false;
            }
        }

        // insert a sieve that optimizes to eliminate composites (this does not appear to be any more efficient than using the hardcoded primes )
        return true;
    }
    
   pub(crate) fn _trial_list(&self, primes: &[u64]) -> bool{
      for i in primes{
        if self.congruence_u64(*i,0){
         return false
        }
      }
      return true
    }

    // weighted strong fermat bases starting from 2^128
    pub(crate) fn weighted_sprp(&self) -> bool {
        const CHECK_LUT: [u8; 5] = [8, 6, 5, 4, 2]; // starting at 2^128    [12, 10, 7, 6, 4, 3, 2, 1];
        let mut check: usize = 1;

        if self.len() < 8 {
            check = CHECK_LUT[self.len() - 3] as usize;
        }

        self.sprp_check(check)
    }

    pub(crate) fn llt(&self, p: u64) -> bool {
        // function will never be called in practical implementation as number-theory does not support the scale of computation needed to use it
        let mut s = Mpz::from_u64(4);

        for _ in 0..(p - 2) {
            s = s.ref_product(&s);
            s.normalize();
            sub_slice(&mut s.limbs[..], &[2]);

            s = s.ref_euclidean(self).1;
        }
        s.normalize();
        if s == Mpz::zero() {
            return true;
        }
        false
    }

    /** Faster than naive evaluation of sophie prime, returns safe prime if true, otherwise None. As this uses is_prime, 
    the accuracy matches that function exactly as the verification of the safe prime itself is deterministic and solely reliant on the accuracy
     of the verification of the sophie prime.
     */
    pub fn is_sophie(&self) -> Option<Self> {
        if self.is_prime() {
            let mut safe = self.shl(1);
            let p = safe.clone();
            let two = Mpz::from_u64(2);
            safe.successor();

            if two.exp_residue(&p, &safe) == Mpz::one() {
                return Some(safe);
            }
        }
        None
    }

    pub(crate) fn jacobi_check_mpz(&self) -> bool {
        // Performs a check of the Jacobi
        let mut witness = 3;
        loop {
            if fast_jacobi_mpz(self, witness) == -1 {
                break;
            }
            witness += 1;
        }

        let witty = Mpz::from_u64(witness);
        self.is_sprp(&witty)
    }

/// Returns an integer in the interval 2^(k-1);2^k  that satisfies the Monier-Rabin bound of passing the Artjuhov-Selfridge test in (1/4)
pub fn psp(k: usize) -> Self{

 let corrector = 1 + (k&1)*2 ;
 
 loop {
  let len = (k-corrector)/2;
  let mut x = Mpz::rand(len);
  if corrector == 3{
   if !x.check_bit(len-2){
    x.flip_bit(len-2);
   }
  }
  if corrector == 1{
    if x.check_bit(len-2){
      x.flip_bit(len-2);
    }
   if x.check_bit(len-3){
      x.flip_bit(len-3);
   }
  }  
  
  let lhs = x.shl(1).ref_addition(&Mpz::one());
  let rhs = x.shl(2).ref_addition(&Mpz::one());
  if lhs.is_prime() & rhs.is_prime(){
   let product = lhs.ref_product(&rhs);
    assert_eq!(product.bit_length(),k as u64);
    
    return product
  }
  }
}
    
  
}

pub(crate) fn fast_jacobi_mpz(x: &Mpz, p: u64) -> i8 {
    let k = x.word_div(p).1;
    k.jacobi(&p)
}

  // detects if the number is in a common pseudoprime form 
pub(crate) fn detect_pseudo_mpz(x: &Mpz) -> bool {
    let mut xminus = x.abs();
    let copy = xminus.clone();
    let eightprod = copy.scale_add(8,1);
    let eightsqrt = eightprod.sqrt().0;
  
    if eightsqrt.sqr() == eightprod{
      return true
    }
  
     let threeprod = copy.scale_add(3,1);
    let threesqrt = threeprod.sqrt().0;
    if threesqrt.sqr() == threeprod{
      return true
    }
    xminus.inv_successor();
    
    
    for i in 1..16 {
        let sq = xminus.word_div(2 * i + 1).0.sqrt().0; // (k*k + k)*(2*i) + k + 1 == x
        let lhs = sq.ref_product(&sq).ref_addition(&sq);
        let lhs2 = lhs.ref_product(&Mpz::from_u64(2 * i + 1));
        if lhs.ref_addition(&lhs2).ref_addition(&Mpz::one()) == copy {
            return true;
        }
    }
    return false;
}
