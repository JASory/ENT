 use std::str;
 use std::str::FromStr;

use crate::arithmetic::sign::Sign;
use crate::arithmetic::sliceops::*;

// converts u64 to string of length 19
fn string_format(x: u64)->String{
    let k = x.to_string();
    let leading = (0..(19-k.len())).map(|_| "0").collect::<String>();
    leading + &k
}

// Naive conversion from Radix-2^64 to Radix-10^19 string, very slow
pub(crate) fn to_string(sign: Sign, digits: Vec<u64>) -> String{

    if digits.len() == 0{
           return "".to_string()
        }
        
    if digits.len() == 1{
       let p = digits[0].to_string();
       if sign == Sign::Negative{
         return "-".to_owned() + &p 
       }
       else{
         return p
       }
    }    
        
       let mut k = vec![];
       let mut x = digits.clone();
       x.push(0u64);
       let xlen = x.len();
       let mut interval = 0usize;
      for i in 0..xlen+1 {
       
         k.push(div_slice(&mut x[..],0x8AC7230489E80000u64));//(xlen-interval)
         
      }
      let mut count=0usize;
      for i in k.iter().rev(){
       if *i > 0u64{
        break;
       }
       count+=1;
      }

      k.truncate(k.len()-count);
        let len = k.len()-1;
        let interim = k[..len].iter().rev().map(|x| string_format(*x)).collect::<Vec<String>>();
        let last = k[len].to_string() + &interim.join("");
        
        if sign == Sign::Negative {
           return "-".to_string() + &last
        }
        return last

} 


pub(crate) fn from_string(string: &str) -> Option<Vec<u64>>{

   let b = string.bytes().all(|c| c.is_ascii_digit());
   
   if b == false{
      return None
   }

     let ln = 19usize;
     let strln = string.len();
     let strem = strln%ln;
     let t = string.as_bytes();
     let mut start = 0usize;
     
     let mut k = vec![];
     
     for i in 0..strln/ln{
         k.push(str::from_utf8(&t[(strln -(i+1)*ln)..(strln-i*ln)]).unwrap()) 
     }
     if strem !=0 { // if there is a remaining slice greater than zero
   if  str::from_utf8(&t[start..strem]).unwrap() != "" { // if the first value is not empty then it pushes it to the vector
          k.push(str::from_utf8(&t[start..strem]).unwrap());
      }  
   }
    let x =  k.iter().rev().map(|x| u64::from_str(x).unwrap()).collect::<Vec<u64>>(); // string converted to decimal vector



       let mut z = vec![0];
        for i in 0..x.len(){
          let mut carry = scale_slice(&mut z[..], 0x8AC7230489E80000);
          if carry > 0 {
            z.push(carry)
          }
          carry = add_slice(&mut z[..], &x[i..(i+1)]) as u64;
          if carry > 0{
            z.push(carry)
          }
        }
        Some(z)
}
