use crate::ntrait::{NumberTheory,MontArith};

pub(crate) fn drbg(mut x: u64) -> u64{
    x ^= x.wrapping_shr(12);
    x ^= x.wrapping_shl(25);
    x ^= x.wrapping_shr(27);
    x.wrapping_mul(0x2545F4914F6CDD1D)
}

 fn poly_eval_32(x: u32, subtrahend: u32, n: u32, npi: u32) -> u32{
   x.mont_prod(x,npi,n).mont_sub(subtrahend,n)
}

pub(crate) fn pollard_brent_32(base: u32,inv:u32,subtrahend: u32, n: u32) -> Option<u32>{
    let m = 4;
    let mut r = 1;
    let mut q = 1;
    let mut g = 1;
    let mut ys = 1;
    let mut y = base;
    let mut x = y;
    
    for _ in 1..19{
      x = y;
      for _ in 0..r{
        y = poly_eval_32(y,subtrahend,n,inv);      
      }
      
      let mut k = 0;
      
      loop{
      
        for i in 0..m{
           if i >= r-k{
           
             break;
           }
         
         y=poly_eval_32(y,subtrahend,n,inv);
         q=q.mont_prod(x.abs_diff(y),inv,n);//machine_prime::mont_prod(q,x.abs_diff(y),n,inv);
         } // end loop

         ys=y;
         g = q.gcd(n);
         k+=m;
         if k >= r || g !=1{
            break;
         }
      }
      
      r<<=1;
      if g != 1{
         break;
      }
      
    }
    
    if g ==n{
       while g==1{
         ys=poly_eval_32(ys,subtrahend,n,inv);
         g=x.abs_diff(ys).gcd(n);
      
       }
    }
    if g !=1 && g !=n && machine_prime::is_prime_wc(g as u64){
       return Some(g);
    }
    None
}

pub(crate) fn poly_factor_32(n: u32) -> u32{
   
      // Start with the polynomial x^2 -1 
      // This works in most cases and is particularly fast to initialise in Montgomery form
   
        let inv = n.inv_2();
        let one = n.n_identity();
        let mut base = one.wrapping_add(one);
        if base > n{
          base = base.wrapping_sub(n);
        }
        
   match pollard_brent_32(base,inv,one,n){
      Some(factor) => return factor,
      None => {
      // if x^2 -1 failed try x^2+1
      // No particular reason except to reuse some values 
        let coef = (n-1).to_mont(n);// machine_prime::to_mont(n-1,n);
        match pollard_brent_32(base,inv,coef,n){
           Some(factor) => return factor,
           None => {
             // Loop that has a roughly 0.5 probability of factoring each iteration
            let mut param = drbg(n.into());
              loop{
                 let  rand_base= (param as u32)%(n-3)+3;
                match pollard_brent_32(rand_base,inv,one,n){
                   Some(factor) => return factor,
                   None => param=drbg(param),
                 }
              }
           }
        }
      }
   }
  }


const fn poly_eval(x: u64, subtrahend: u64, n: u64, npi: u64) -> u64{
     machine_prime::mont_sub(machine_prime::mont_prod(x,x,n,npi),subtrahend,n)
}

fn pollard_brent(base: u64,inv:u64,subtrahend: u64, n: u64) -> Option<u64>{
    let m = 128;
    let mut r = 1;
    let mut q = 1;
    let mut g = 1;
    let mut ys = 1;
    let mut y = base;
    let mut x = y;
    
    for cycle in 1..17{    
      x = y;
      for _ in 0..r{
        y = poly_eval(y,subtrahend,n,inv);      
      }
      
      let mut k = 0;
      
      loop{
      
        for i in 0..(m*cycle){
           if i >= r-k{
             break;
           }
         
         y=poly_eval(y,subtrahend,n,inv);
         q=machine_prime::mont_prod(q,x.abs_diff(y),n,inv);
         } // end loop

         ys=y;
         g = q.gcd(n);
         k+=m;
         if k >= r || g !=1{
            break;
         }
      }
      
      r<<=1;
      if g != 1{
         break;
      }
      
    }
    
    if g ==n{
       while g==1{
         ys=poly_eval(ys,subtrahend,n,inv);
         g=x.abs_diff(ys).gcd(n);
      
       }
    }
    if g !=1 && g !=n && machine_prime::is_prime_wc(g){
       return Some(g);
    }
    None
}

pub(crate) fn poly_factor(n: u64) -> u64{
   
      // Start with the polynomial x^2 -1 
      // This works in most cases and is particularly fast to initialise in Montgomery form
   let inv = machine_prime::mul_inv2(n);
   let one = machine_prime::one_mont(n);
   let base = machine_prime::two_mont(one,n);
   
   match pollard_brent(base,inv,one,n){
      Some(factor) => return factor,
      None => {
      // if x^2 -1 failed try x^2+1
      // No particular reason except to reuse some values 
        let coef = machine_prime::to_mont(n-1,n);
        match pollard_brent(base,inv,coef,n){
           Some(factor) => return factor,
           None => {
             // Loop that has a roughly 0.5 probability of factoring each iteration
            let mut param = drbg(n);
              loop{
                 let  rand_base= param%(n-3)+3;
                match pollard_brent(rand_base,inv,one,n){
                   Some(factor) => return factor,
                   None => param=drbg(param),
                 }
              }
           }
        }
      }
   }
  }


const fn poly_eval_128(x: u128, subtrahend: u128, n: u128, npi: u128) -> u128{
     machine_prime::mont_sub_128(machine_prime::mont_sqr_128(x,n,npi),subtrahend,n)
}

fn pollard_brent_128(base: u128,inv:u128,subtrahend: u128, n: u128) -> Option<u128>{
    let m = 512;
    let mut r = 1;
    let mut q = 1;
    let mut g = 1;
    let mut ys = 1;
    let mut y = base;
    let mut x = y;
    
    for cycle in 1..34{    
      x = y;
      for _ in 0..r{
        y = poly_eval_128(y,subtrahend,n,inv);      
      }
      
      let mut k = 0;
      
      loop{
      
        for i in 0..(m*cycle){
           if i >= r-k{
             break;
           }
         
         y=poly_eval_128(y,subtrahend,n,inv);
         q=machine_prime::mont_prod_128(q,x.abs_diff(y),n,inv);
         } // end loop

         ys=y;
         g = q.gcd(n);
         k+=m;
         if k >= r || g !=1{
            break;
         }
      }
      
      r<<=1;
      if g != 1{
         break;
      }
      
    }
    
    if g ==n{
       while g==1{
         ys=poly_eval_128(ys,subtrahend,n,inv);
         g=x.abs_diff(ys).gcd(n);
      
       }
    }
    if g !=1 && g !=n && machine_prime::is_prime_wc_128(g){
       return Some(g);
    }
    None
}




 pub(crate) fn poly_factor_128(n: u128) -> u128{
   
      // Start with the polynomial x^2 -1 
      // This works in most cases and is particularly fast to initialise in Montgomery form
   let inv = machine_prime::mul_inv2_128(n);
   let one = machine_prime::one_mont_128(n);
   let base = machine_prime::two_mont_128(one,n);
   
   match pollard_brent_128(base,inv,one,n){
      Some(factor) => return factor,
      None => {
      // if x^2 -1 failed try x^2+1
      // No particular reason except to reuse some values 
        let coef = machine_prime::to_mont_128(n-1,n);
        match pollard_brent_128(base,inv,coef,n){
           Some(factor) => return factor,
           None => {
             // Loop that has a roughly 0.5 probability of factoring each iteration
            let mut param = drbg(n as u64);
              loop{
                 let  rand_base= (param as u128)%(n-3)+3;
                match pollard_brent_128(rand_base,inv,one,n){
                   Some(factor) => return factor,
                   None => param=drbg(param),
                 }
              }
           }
        }
      }
   }
  }

