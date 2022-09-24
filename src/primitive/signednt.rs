use crate::traits::NumberTheory;

macro_rules! signednt(
  ($($t:ty; $s:ty),* $(,)*) => {$(
impl NumberTheory for $t{

       fn rng() -> $t {
         <$s>::rng() as $t
       }

       fn euclidean_div(&self, other: &Self) -> (Self, Self) {
           (*self/ *other, *self % *other)
       }

       fn is_sprp(&self, base: &Self) -> bool {
        if base < &0{
         if base > self {
          return (self.abs() as $s).is_sprp(&((self.abs()+(base % self.abs())) as $s))
         }
          return (self.abs() as $s).is_sprp(&((self.abs() + *base) as $s))
        }
           (self.abs() as $s).is_sprp(&(*base as $s))
       }

       fn is_prime(&self) -> bool {
          (self.abs() as $s).is_prime()
       }

      fn prime_proof(&self) ->(bool, Vec<Self>) {
         let (flag, certificate) = (self.abs() as $s).prime_proof();
         (flag,certificate.iter()
          .map(|x| *x as $t)
          .collect::<Vec<$t>>())
       }

       fn prime_list(&self, sup: &Self) -> Vec<Self> {

        let inf = std::cmp::min(*self, *sup);
        let mut hi = std::cmp::max(*self,*sup);

        if hi < Self::MAX{
           hi+=1;
        }

        let mut primevector = vec![];

        for i in inf..hi{
          if i.is_prime(){
           primevector.push(i)
          }
        }
        primevector
    }

   fn nth_prime(&self) -> Option<Self> {
      let mut count : Self = 0;
      let mut start : Self =  0;
      loop {
         start+=1;
        if start == Self::MAX{
          return None
        }
        if start.is_prime(){
          count+=1;
        }
        if count == *self{
          return Some(start)
        }
      }
    }

    fn pi(&self) -> Self{
      let mut count : Self =  0;
      let zero : Self = 0;
      for i in zero..*self{
       if i.is_prime(){
         count+=1;
       }
      }
      count
    }

    fn prime_gen(x: u32) -> Self{
      if x > Self::BITS-1{
        panic!("Does not fit into datatype")
      }
      <$s>::prime_gen(x) as $t
    }

       fn factor(&self) -> Vec<Self> {
         (self.abs() as $s)
          .factor().iter()
          .map(|x| *x as $t)
          .collect::<Vec<$t>>()
       }

       fn sqrt(&self) -> (Self, Self) {
         if *self < 0 {
           return ((self.abs() as $s).sqrt().0 as $t,1)
         }
         ((*self as $s).sqrt().0 as $t,0)
       }

       fn nth_root(&self, n: &Self) -> (Self, Self) {
         if *n < 0{return (0,0)}

         if  *n == 1{
           return (*self,0)
         }

         if *self < 0 && n.abs()%2 == 0{
           return ((self.abs() as $s).nth_root(&(n.abs() as $s)).0 as $t,1)
         }
         return ((self.abs() as $s).nth_root(&(n.abs() as $s)).0 as $t,0)


       }

       fn radical(&self) -> Self {
         (self.abs() as $s).radical() as $t
       }

       fn k_free(&self, k: &Self) -> bool {
         (self.abs() as $s).k_free(&(k.abs() as $s))
       }

       fn euclid_gcd(&self, other: &Self) -> Self {
        (self.abs() as $s).euclid_gcd(&(other.abs() as $s)) as $t
       }


       fn extended_gcd(&self, other: &Self)->(Self,Self,Self){
         let mut gcd : Self =*self;
         let mut new_r : Self =*other;
         let mut bezout_1 : Self =1;
         let mut new_s : Self =0;
         let mut bezout_2 : Self = 0;
         let mut new_t: Self = 1;

    while new_r !=0 {
    let quotient =gcd/new_r;
    let mut temp : Self =new_r;
    new_r=gcd-quotient*temp;
    gcd=temp;

    temp=new_s;
    new_s=bezout_1-quotient*temp;
    bezout_1=temp;

    temp=new_t;
    new_t=bezout_2-quotient*temp;
    bezout_2=temp;

    }
    (gcd,bezout_1,bezout_2)
}

       fn lcm(&self, other: &Self) -> Self{
         (self.abs() as $s).lcm(&(other.abs() as $s)) as $t
       }

       fn checked_lcm(&self, other: &Self) -> Option<Self> {
        (self.abs() as $s).checked_lcm(&(other.abs() as $s)).map(|y| y as $t)
       }

       fn euler_totient(&self) -> Self {
         (self.abs() as $s).euler_totient() as $t
       }

       fn jordan_totient(&self, k: &Self) -> Option<Self> {
          (self.abs() as $s).jordan_totient(&(k.abs() as $s)).map(|y| y as $t)
       }

       fn dedekind_psi(&self, k: &Self) -> Option<Self> {
          (self.abs() as $s).dedekind_psi(&(k.abs() as $s)).map(|y| y as $t)
       }

       fn quadratic_residue(&self, n: &Self) -> Self {
        (self.abs() as $s).quadratic_residue(&(n.abs() as $s)) as $t
       }

       fn product_residue(&self, other: &Self, n: &Self) -> Self {
        let mut a = self.clone();
        let mut b = other.clone();
        let modulo = n.abs();

        if a < 0 {
            a += modulo;
        }
        if b < 0 {
            b += modulo;
        }
        (a as $s).product_residue(&(b as $s), &(modulo as $s)) as $t
    }

    fn exp_residue(&self, pow: &Self, n: &Self) -> Self {
        let mut a = self.clone();
        if a < 0 {
            a = n.abs() + self
        }
        if pow < &0{
          let (gcd, inv,_) = (a as $s).extended_gcd(&(n.abs() as $s));
          if gcd != 1{
           panic!("No multiplicative inverse");
          }
          return inv.exp_residue(&(pow.abs() as $s), &(n.abs() as $s)) as $t
        }
        (a as $s).exp_residue(&(pow.abs() as $s), &(n.abs() as $s)) as $t
    }

    fn checked_exp_residue(&self, pow: &Self, n: &Self) -> Option<Self> {
        let mut a = self.clone();

        if a < 0 {
            a = n.abs() + self
        }

        if pow < &0{
          let (gcd, inv,_) = (a as $s).extended_gcd(&(n.abs() as $s));
          if gcd != 1{
           return None
          }
          return inv.checked_exp_residue(&(pow.abs() as $s), &(n.abs() as $s)).map(|y| y as $t) //as $t
        }
        (a as $s).checked_exp_residue(&(pow.abs() as $s), &(n.abs() as $s)).map(|y| y as $t)
    }

    fn legendre(&self, p: &Self) -> i8 {
        let k = self.exp_residue(&((p.abs() - 1) >> 1), &p.abs());
        if k == 1 {
            return 1;
        };
        if k == p.abs() - 1 {
            return -1;
        };
        return 0;
    }

    fn checked_legendre(&self, p: &Self) -> Option<i8> {
        if p.abs() == 2 || p.is_prime() == false {
            return None;
        }
        Some(self.legendre(&p))
    }

    fn liouville(&self) -> i8{
      (*self as $s).liouville()
    }
    
    fn derivative(&self) -> Option<Self> {
      (self.abs() as $s).derivative().map(|y| y as $t)
    }

    fn mangoldt(&self) -> f64 {
      (self.abs() as $s).mangoldt()
    }
    
    fn mobius(&self) -> i8 {
      (self.abs() as $s).mobius()
    }

    fn jacobi(&self, k: &Self) -> i8 {

        let mut n = *self;
        let mut p = *k;
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

    fn checked_jacobi(&self, k: &Self) -> Option<i8> {
        if k > &0 && *k % 2 == 1 {
            return Some(self.jacobi(k));
        }
        return None;
    }

    fn smooth(&self) -> Self {
        let k = self.factor();
        k[k.len() - 2]
    }

    fn is_smooth(&self, b: &Self) -> bool {
        self.smooth() <= b.abs()
    }
  }
    )*}
);

signednt!(i8;u8, i16;u16, i32;u32, i64;u64, isize;usize, i128;u128);
