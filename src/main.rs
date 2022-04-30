use number_theory::NumberTheory;
use number_theory::Mpz;
use number_theory::Sign;

const PRIMES : [u32;128] = [
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
   73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
   179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281,
   283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409,
   419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541,
   547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659,
   661, 673, 677, 683, 691, 701, 709, 719, 727,
  
   ];
 

fn sprp_32(p: u32, base: u32)->bool{// checks if base^p = 1 mod p  or base^(d*2^n)= -1 for some n  
     let zeroes = (p-1).trailing_zeros() as u32; // Breaks number down to p= d*2^n -1
     let d = (p-1)/ (1<<zeroes);
     let mut x = base.mod_pow(&d,&p); // base^d mod p
     if x == 1u32 || x==p-1{   // checks if base^p = 1 mod p  or base^(d*2^n)= -1
       return true
       }
    for _ in 0..zeroes-1{// checks for all d*2^zeroes. One is subtracted since d*2^n was already checked above
     x = x.quadratic_residue(&p);
     if x == p-1 {       // if any d*2^zeroes = p-1  then it passes
       return true
     }
    }
    return false        // otherwise it fails
 }
 
 fn sprp_64(p: u64, base: u64)->bool{// checks if base^p = 1 mod p  or base^(d*2^n)= -1 for some n  
     let zeroes = (p-1).trailing_zeros() as u64; // Breaks number down to p= d*2^n -1
     let d = (p-1)/ (1<<zeroes);
     let mut x = base.mod_pow(&d,&p); // base^d mod p
     if x == 1u64 || x==p-1{   // checks if base^p = 1 mod p  or base^(d*2^n)= -1
       return true
       }
    for _ in 0..zeroes-1{// checks for all d*2^zeroes. One is subtracted since d*2^n was already checked above
     x = x.quadratic_residue(&p);
     if x == p-1 {       // if any d*2^zeroes = p-1  then it passes
       return true
     }
    }
    return false        // otherwise it fails
 }
 
 fn trial_prime(x: u32, base: u32) -> bool{
    for i in PRIMES.iter(){
      if x == *i as u32 {return true}
      if x % *i as u32 == 0 {return false}
    }
     sprp_32(x,base)
    }
    
 fn trial_prime64(x: u64, base: u64) -> bool{
    for i in PRIMES.iter(){
      if x == *i as u64 {return true}
      if x % *i as u64 == 0 {return false}
    }
     sprp_64(x,base)
 }   
 
 fn sqrt_precision(x: Mpz, places: u64, linelen: usize) -> String {
      let lead = x.sqrt().to_string();
      let time = std::time::Instant::now();
       let dec = Mpz::from_u64(5).pow(places*2).shl((places*2) as usize);
       let stop = time.elapsed();
      // println!("{:?}", stop);
      let tmp = x.ref_product(&dec).sqrt().to_string(); 
      let len = tmp.len();
      let trailing = tmp.as_bytes();

           let mut result : Vec<u8> = vec![];
           for i in lead.as_bytes(){
             result.push(*i)
           }
           result.push(46); // inserts decimal point  46  10
           let len = trailing.len();
           for i in 0..len{//.as_bytes()[lead.len()..].iter(){
            if i%linelen == 0 && i > 0{
              result.push(10)
            }
             result.push(trailing[i])
           }
           
           String::from_utf8(result).unwrap()
 }

fn main(){/*
  let mut bound = 2000u64;
  let time = std::time::Instant::now();
  for k in 1..241{
 'outer : for j in 2..65535{
     let mut count = 0u64;
  'inner : for i in ((1u64<<32) + 20_000_000)..((1u64<<32) + 40_000_000){
    if i.is_prime() != trial_prime64(i,j) {
    //println!("{}",i);
      count+=1;
      if count >= bound{
        break 'inner ;
      } 
    }
  }
   if bound > count{
     println!("{} {}",j,count);
     bound = count;
     count = 0u64;
   }
   
  }
  let stop = time.elapsed();
  println!("{:?} primes counted in {:?}",bound,stop)
  }
  /*
    let mut point = 0u32;
  for i in 4194304..u32::MAX{
     if (((i-1) as f64).ln()*30f64) as usize != (((i) as f64).ln()*30f64) as usize{
       println!("{},",i)
     }
  }
  */
}
*/
/*
for k in 1..215{
 let inf = ((1u64<<32) + k*20_000_000);
 let sup = ((1u64<<32) + (k+1)*20_000_000);

 'outer : for j in 2..65535{

   let mut count = 0u64;
  'inner : for i in inf..sup{
    if i.is_prime() != trial_prime64(i,j) {
       count +=1;
      break 'inner;
    }
  }
   if count == 0{
     println!("{}-base",j);
      break 'outer;
    }
  }
  }

}
*/
 /*
 for i in (1u64<<32)+40_000_000..(1u64<<32)+60_000_000{
    if i.is_prime() != trial_prime64(i,20){
      println!("Error!")
    }
 }
 }
 */
 /*
 let time = std::time::Instant::now();
 let mut count = 0u64;
 for i in 1u64<<32..(1u64<<33){
  if i.is_prime(){
  //println!("Executed");
   count+=1;
  }
 }
 let elapsed = time.elapsed();
 println!("{}primes counted in {:?}",count,  elapsed);
 */
 const PRECISION : u64 = 3_000_000;
 let two = Mpz::from_u64(2);
 let digits = Mpz::from_u64(5).pow(PRECISION).shl(PRECISION as usize);
 let total = two.ref_product(&digits);
 let k = total.nth_root(3);
 let z = k.pow(3);
 //let p = k.ref_product(&k);
 //println!("{}",total.to_string());
 println!("{:?}",z.to_string())
 }
