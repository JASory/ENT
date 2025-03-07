use crate::ntrait::NumberTheory;

/// A primality certificate, currently a modification of the Pratt Certificate
/// Instead of checking if a^(n-1) mod n we use the stronger form of the Fermat test
/// which permits proving that Carmichael numbers are composite  
#[derive(Clone)]
pub struct Certificate<T: NumberTheory>{
   pub(crate) n: T,
   pub(crate) witness: T,
   pub(crate) fctr: Vec<T>,
}

impl<T: NumberTheory> Certificate<T>{

  pub(crate) fn new(n: T, witness: T, fctr: Vec<T>) -> Self{
     Self{n,witness,fctr}
  } 
  
 
    /// Check the certificate 
   pub fn check(&self) -> bool{
     if self.fctr.is_empty(){
        return false;
     }

     if !self.n.clone().strong_fermat(self.witness.clone()){
         return false;
     }

     for i in self.fctr.iter(){
     if self.witness.exp_residue(i.clone(), self.n.clone()).is_unit() {
         return false
        }
     }
   return true;
   }
   
}

#[derive(Clone,Default)]
pub struct Factorization<T>{
   pub(crate) base: Vec<T>,
   pub(crate) power: Vec<u64>,
}

impl <T: NumberTheory> Factorization<T>{

   pub fn new() -> Self{
      Self{base: vec![], power: vec![]}
   }
   
   pub fn from_components(base: Vec<T>,power: Vec<u64>) -> Self{
        Self{base,power}
   }
   
   pub fn factor_iter(&self) -> std::slice::Iter<T>{
      self.base.iter()
   }
   
   pub fn power_iter(&self) -> std::slice::Iter<u64>{
     self.power.iter()
   }
   
   pub fn add_factor(&mut self, fctr: T){
      self.base.push(fctr);
   }
   
   pub fn add_power(&mut self, power: u64){
      self.power.push(power);
   }
   
  pub fn pair_iter(&self) -> std::iter::Zip<std::slice::Iter<T>,std::slice::Iter<u64>>{
      self.base.iter().zip(self.power.iter())
  } 
  
  pub fn max(&self) -> T{
     self.factor_iter().max().unwrap().clone() 
  }
  
  pub fn k_free(&self, k: u64) -> bool{
      for i in self.power_iter(){
        if *i >= k{
           return false;
        }
      }
      return true;
  }
  
  pub fn prime_omega(&self) -> u64{
     self.power_iter().sum::<u64>()
  }
  
  }
  
 impl<T: NumberTheory> std::fmt::Display for Factorization<T>{
 
   fn fmt(&self,f: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error>{

      let mut output = String::new();
      if self.power[0] != 1{
        output+=&(self.base[0].to_string()+"^"+&self.power[0].to_string());
      }
      else{
        output+=&(self.base[0].to_string());
      }
      if self.base.len() > 1{
      for i in 1..self.base.len(){
         if self.power[i] != 1{
         let pair = self.base[i].to_string()+"^"+&self.power[i].to_string();
         output+=&(" * ".to_owned() + &pair);
         }
         else{
         let pair = self.base[i].to_string();
         output+= &(" * ".to_owned()+&pair);
         }
      }
      }
      write!(f,"{}",output)
  }
}


