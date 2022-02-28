use crate::arithmetic::sliceops::*;
use crate::arithmetic::inlineops::*;
use crate::arithmetic::sign::Sign;
use std::cmp::Ordering;



pub(crate) fn sub_sign( mut a: &[u64], mut b: &[u64]) -> (Sign,Vec<u64>) {
    if let Some(&0) = a.last() {
        a = &a[..a.iter().rposition(|&x| x != 0).map_or(0, |i| i + 1)];
    }
    if let Some(&0) = b.last() {
        b =  &b[..b.iter().rposition(|&x| x != 0).map_or(0, |i| i + 1)];
    }

    match cmp_slice(a, b) {
        Ordering::Greater => {
           let mut a = a.to_vec();
            sub_slice( &mut a[..], b);   
            (Sign::Positive,a)  
        }
        Ordering::Less => {
          let mut b = b.to_vec();
            sub_slice( &mut b[..], a);
            (Sign::Negative,b)
        }
        Ordering::Equal => (Sign::Positive,vec![0]),
    }
}

pub(crate) fn mac_scalar(x: &[u64], y: u64,acc: &mut[u64]) {
    if y == 0 {
        return;
    }

    let mut carry = 0;
    let (a_lo, a_hi) = acc.split_at_mut(x.len());

    for (a, &x) in a_lo.iter_mut().zip(x) {
        *a = mul_acc(*a, x, y, &mut carry);
    }
   
   let (carry_hi,carry_lo) =  split(carry);

    let final_carry = if carry_hi == 0 {
        add_slice(a_hi, &[carry_lo])
    } else {
        add_slice(a_hi, &[carry_hi, carry_lo])
    };

}
  // scale and subtract y from x
fn sub_mul(x: &mut [u64], y: &[u64], z: u64)->u64{
    let mut carry = u64::MAX;
    
    for (i,j) in x.iter_mut().zip(y){
        let fused = fuse(u64::MAX, *i);
        let result = fused - u64::MAX as u128 + carry as u128 - *j as u128 * z as u128;

        let (new_x, new_carry) = split(result);
        carry = new_carry;
        *i = new_x;
    }
    u64::MAX-carry
}


pub(crate) fn mul_slice(mut b: &[u64], mut c: &[u64],mut acc: &mut[u64]){

if let Some(&0) = b.first() {
        if let Some(nz) = b.iter().position(|&d| d != 0) {
            b = &b[nz..];
            acc = &mut acc[nz..];
        } else {
            return;
        }
    }
    if let Some(&0) = c.first() {
        if let Some(nz) = c.iter().position(|&d| d != 0) {
            c = &c[nz..];
            acc = &mut acc[nz..];
        } else {
            return;
        }
    }

    let acc = acc;
    let (x, y) = if b.len() < c.len() { (b, c) } else { (c, b) };
    
   if x.len() <= 32 {
    for (i, xi) in x.iter().enumerate() {
            mac_scalar(y, *xi,&mut acc[i..]);
        }
    }
    
    else {
    
        let b = x.len() / 2;
        let (x0, x1) = x.split_at(b);
        let (y0, y1) = y.split_at(b);

        
        let len = x1.len() + y1.len() + 1;
        let mut p =  vec![0; len];


        mul_slice(&x1[..], &y1[..],&mut p[..]);

       
        remove_lead_zeros(&mut p);

        add_slice(&mut acc[b..], &p);
        add_slice(&mut acc[b * 2..], &p);


        p.truncate(0);
        p.resize(len, 0);


        mul_slice(x0, y0,&mut p);
        remove_lead_zeros(&mut p);

        add_slice(acc, &p);
        add_slice(&mut acc[b..], &p);

        
        let (mut j0_sign, j0) = sub_sign(x1, x0);
        let (j1_sign, j1) = sub_sign(y1, y0);

        match j0_sign.mul(&j1_sign) {
         Sign::Positive => {
                p.truncate(0);
                p.resize(len, 0);

                mul_slice(&j0[..], &j1[..],&mut p[..]);
                remove_lead_zeros(&mut p);

                sub_slice(&mut acc[b..], &p);
            }
            Sign::Negative => {
                mul_slice(&j0, &j1,&mut acc[b..]);
            }
        }
    
    
    }    
}
 
 
pub(crate) fn euclidean_slice(a: &mut Vec<u64>, b: &[u64])->(Vec<u64>,Vec<u64>){
 
  let mut a0 = 0;


    let b0 = *b.last().unwrap();
    let b1 = b[b.len() - 2];

    let quo_len = a.len() - b.len() + 1;
    let mut quo =  vec![0; quo_len];

     for j in (0..quo_len).rev() {
       
        let a1 = *a.last().unwrap();
        let a2 = a[a.len() - 2];
        
         let mut q0 = divide3by2(a0,a1,a2,b0,b1);

        let mut borrow = sub_mul(&mut a[j..], b, q0);

        if borrow > a0 {
         
            q0 -= 1;
            borrow -= add_slice(&mut a[j..], b) as u64;
        }


        quo[j] = q0;

        a0 = a.pop().unwrap();
    }

    a.push(a0);
    let mut remainder = a.to_vec();
    let mut quotient = quo.to_vec();
    remove_lead_zeros(&mut remainder);

     remove_lead_zeros(&mut quotient);
     if remainder.len() == 0{
       remainder.push(0u64)
     }
     
     if quotient.len() == 0{
       
       quotient.push(0u64)
     }
    (quotient,remainder )
}
