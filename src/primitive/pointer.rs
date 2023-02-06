use crate::traits::NumberTheory;

impl NumberTheory for usize {
    fn rng() -> Self {
        u64::rng() as usize
    }
    
    fn residue(&self, ring: &Self) -> Self {
      (*self as u64).residue(&(*ring as u64)) as usize
    }
    fn euclidean_div(&self, other: &Self) -> (Self, Self) {
        let (quo, rem) = (*self as u64).euclidean_div(&(*other as u64));
        (quo as usize, rem as usize)
    }

    fn is_sprp(&self, other: &Self) -> bool {
        (*self as u64).is_sprp(&(*other as u64))
    }

    fn is_prime(&self) -> bool {
        (*self as u64).is_prime()
    }

    fn prime_proof(&self) -> (bool, Vec<Self>) {
        let (flag, cert) = (*self as u64).prime_proof();
        (
            flag,
            cert.iter().map(|x| *x as usize).collect::<Vec<usize>>(),
        )
    }

    fn prime_list(&self, sup: &Self) -> Vec<Self> {
        let list = (*self as u64).prime_list(&(*sup as u64));
        list.iter().map(|x| *x as usize).collect::<Vec<usize>>()
    }

    fn nth_prime(&self) -> Option<Self> {
        (*self as u64).nth_prime().map(|y| y as usize)
    }

    fn pi(&self) -> Self {
        (*self as u64).pi() as usize
    }

    fn prime_gen(x: u32) -> Option<usize> {
        u64::prime_gen(x).map(|z| z as usize)
    }

    fn factor(&self) -> Vec<Self> {
        let list = (*self as u64).factor();
        list.iter().map(|x| *x as usize).collect::<Vec<usize>>()
    }

    fn checked_factor(&self) -> Option<Vec<Self>>{
      if *self < 2{
       return None
      }
      Some(self.factor())
    }


    fn sqrt(&self) -> (Self, Self) {
        ((*self as u64).sqrt().0 as usize, 0)
    }

    fn nth_root(&self, n: &Self) -> (Self, Self) {
        ((*self as u64).nth_root(&(*n as u64)).0 as usize, 0)
    }

    fn euclid_gcd(&self, other: &Self) -> Self {
        (*self as u64).euclid_gcd(&(*other as u64)) as usize
    }

    fn extended_gcd(&self, other: &Self) -> (Self, Self, Self) {
        let (gcd, s_inv, o_inv) = (*self as u64).extended_gcd(&(*other as u64));
        (gcd as usize, s_inv as usize, o_inv as usize)
    }

    fn lcm(&self, other: &Self) -> Self {
        (*self as u64).lcm(&(*other as u64)) as usize
    }

    fn checked_lcm(&self, other: &Self) -> Option<Self> {
        (*self as u64)
            .checked_lcm(&(*other as u64))
            .map(|y| y as usize)
    }

    fn euler_totient(&self) -> Self {
        (*self as u64).euler_totient() as usize
    }

    fn jordan_totient(&self, k: &Self) -> Option<Self> {
        (*self as u64)
            .jordan_totient(&(*k as u64))
            .map(|y| y as usize)
    }
    
    fn carmichael_totient(&self) -> Option<Self>{
      (*self as u64).carmichael_totient().map(|y| y as usize)
    }

    fn dedekind_psi(&self, k: &Self) -> Option<Self> {
        (*self as u64)
            .dedekind_psi(&(*k as u64))
            .map(|y| y as usize)
    }

    fn product_residue(&self, other: &Self, n: &Self) -> Self {
        (*self as u64).product_residue(&(*other as u64), &(*n as u64)) as usize
    }

    fn checked_product_residue(&self, other: &Self, n: &Self) -> Option<Self> {
        (*self as u64).checked_product_residue(&(*other as u64), &(*n as u64)).map(|y| y as usize)
    }
    
    fn quadratic_residue(&self, n: &Self) -> Self {
        (*self as u64).quadratic_residue(&(*n as u64)) as usize
    }

    fn checked_quadratic_residue(&self, n: &Self) -> Option<Self> {
        (*self as u64).checked_quadratic_residue(&(*n as u64)).map(|y| y as usize)
    }
    
    fn exp_residue(&self, pow: &Self, n: &Self) -> Self {
        (*self as u64).exp_residue(&(*pow as u64), &(*n as u64)) as usize
    }

    fn checked_exp_residue(&self, pow: &Self, n: &Self) -> Option<Self> {
        (*self as u64)
            .checked_exp_residue(&(*pow as u64), &(*n as u64))
            .map(|y| y as usize)
    }

    fn k_free(&self, k: &Self) -> bool {
        (*self as u64).k_free(&(*k as u64))
    }


    fn radical(&self) -> Option<Self>{
        (*self as u64).radical().map(|y| y as usize)
    }

    fn smooth(&self) -> Self {
        (*self as u64).smooth() as usize
    }

    fn is_smooth(&self, b: &Self) -> bool {
        (*self as u64).is_smooth(&(*b as u64))
    }

    fn legendre(&self, p: &Self) -> i8 {
        (*self as u64).legendre(&(*p as u64))
    }

    fn liouville(&self) -> i8 {
        (*self as u64).liouville()
    }
    
    fn derivative(&self) -> Option<Self>{
       (*self as u64).derivative().map(|y| y as usize)
    }

    fn checked_legendre(&self, p: &Self) -> Option<i8> {
        (*self as u64).checked_legendre(&(*p as u64))
    }

    fn mangoldt(&self) -> f64 {
        (*self as u64).mangoldt()
    }
    
    fn mobius(&self) -> i8 {
      (*self as u64).mobius()
    }


    fn jacobi(&self, p: &Self) -> i8 {
        (*self as u64).jacobi(&(*p as u64))
    }

    fn checked_jacobi(&self, p: &Self) -> Option<i8> {
        (*self as u64).checked_jacobi(&(*p as u64))
    }

}
