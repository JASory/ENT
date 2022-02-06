/*

  Primality 

*/
  use crate::Mpz;
  use crate::traits::NumberTheory;
  use crate::primes::MERSENNE_LIST;
  use crate::primes::PRIMELIST;
  
  use crate::arithmetic::sliceops::sub_slice;
  use crate::arithmetic::inlineops::rng_64;
  use crate::arithmetic::sliceops::mod_slice;



 impl Mpz {
 
 
  pub fn sprp(&self, base: Self)->bool{
      let mut p_minus = self.clone();
      let one = Mpz::one();

      sub_slice(&mut p_minus.limbs[..],&one.limbs[..]); //subtract one this will always succeed

      let zeroes = p_minus.trailing_zeros() as usize;

      let d = p_minus.shr(zeroes);  
      let mut x = base.u_mod_pow(&d, self);
      
      if x == Mpz::one() || x == p_minus {
         return true
      }
      
      for i in 0..zeroes -1{

        x = x.u_quadratic_residue(&self);

        if x == p_minus{
          return true
        }
      }
      return false
  }
      
  pub fn sprp_check(&self, steps: usize) -> bool{
      if self.len() < 2 {return self.to_u64().unwrap().is_prime()}  // if fits into u64, reduce to 64-bit check 
    
      if self.is_fermat(){return false}
      
      for i in 0..steps{
     let z = Mpz::rand(self.len(),rng_64).ref_euclidean(&self).1;

     if self.sprp(z)==false{ return false} //
   }
   
   return true
  } 
  
    // weighted to maintain at most 2^-64 probability of failure, slower than most implementations for small numbers but faster for larger. Values greater than 2^512 receive only two checks, a strong-base 2 and a random-base check. This is due to the fact that the density of pseudoprimes rapidly declines

  pub fn probable_prime(&self) -> bool{
      const CHECK_LUT : [u8;10]  = [12,11,9,6,5,5,4,3,2,1];
      const DIV_BOUND : [u16;10] = [380,500,800,1000,1200,1400,1600,1800,2000,2048]; // 8 300, 16 500  32 1000  64 2048
      let mut check = 1;
      let mut supremum = 2048;
      if self.len() < 2usize {
         return self.to_u64().unwrap().is_prime()
      }
      
     if self.len() < 12{
         check = CHECK_LUT[self.len()-2];
         supremum = DIV_BOUND[self.len()-2];
     }

   let two = Mpz::from_u64(2);

   if self.is_even(){return false}
   if self.is_fermat(){return false}
   

   
    match self.is_mersenne() {
      Some(x) => {if MERSENNE_LIST.contains(&(x as u32)){return true}
                   else if x < 57885161 {return false}
                   else{return self.llt(x)} }
      None    => (),
    }
   
    let rem = mod_slice(&self.limbs[..],16294579238595022365);
    
    for i in PRIMELIST[1..16usize].iter(){
        if rem%*i as u64 == 0{
          return false
        }
    }
    
  // println!("{}", supremum);
   
   for i in PRIMELIST[17..2048u64 as usize].iter(){ // 295 14.074467005s   295 13.65539096s  


     if self.congruence_u64(*i as u64,0){
       return false
     }
   }
   
   //println!("{}",check);
   if self.sprp(two)==false{return false}
 
  let z = self.sprp_check(check as usize +2);
     
   return z
  
  }
  
  pub(crate) fn llt(&self, p: u64) -> bool{// function will never be called in practical implementation
  	let mut s = Mpz::from_u64(4);

  	
  	for i in 0..(p-2){
  	  s = s.ref_product(&s);
  	  s.normalize();
  	  sub_slice(&mut s.limbs[..], &[2]);
  	  
  	  s = s.ref_euclidean(&self).1;
  	}
  	s.normalize();
  	if s == Mpz::zero(){return true}
  	return false
  }
 
 
  pub fn is_sophie(&self) -> Option<Self> {
  	
  	if self.is_prime(){
  	  let mut safe = self.shl(1);
  	  let p= safe.clone();
  	  let two = Mpz::from_u64(2);
  	  safe.successor();
  	 
  	 if two.mod_pow(&p,&safe) == Mpz::one() {
  	    return Some(safe)
  	 }
  	 
  	}
  	return None
  }
 
 }
