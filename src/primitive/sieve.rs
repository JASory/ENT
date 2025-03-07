
pub(crate) fn sieve(sup: u64, collect: bool) -> (Vec<u64>,u64){ 
   let isqrt = (sup as f64).sqrt() as u64;
   let segment_size = std::cmp::max(isqrt,32768);
   let mut count = 0u64;
   let mut values = vec![];
   let mut sieve = Vec::<bool>::with_capacity(segment_size as usize);
   let mut is_prime = vec![true;isqrt as usize+1];
   let mut primes = vec![];
   let mut multiples = vec![];
   let mut low = 0u64;
   
   let mut i : u64 = 3;
   let mut n : u64 = 3;
   let mut s : u64 = 3; 
   
   loop { // full loop
     
    
     if low > sup{
       break;
     }
     
     
     sieve.truncate(0);
     sieve.resize(segment_size as usize,true); // allocate with true values
     
     let mut high : u64 = low + segment_size -1;
     high = std::cmp::min(high,sup);
     
     let inner_sup = (high as f64).sqrt() as u64;
    
     loop{ // Generate sieving primes
     
       if i*i > high{
          break;
       }
       
       if is_prime[i as usize]{
                      let mut j = i*i;
         loop {

            if j > isqrt{
               break;
            }
            is_prime[j as usize] = false;
                        j+=i;
         }
         
       }
              i+=2;
     } // End prime generation
    
     loop{// prime initialisation
       if s*s > high{
         break;
       }
       
       if is_prime[s as usize]{
          primes.push(s);
          multiples.push(s*s -low);
       }
              s+=2;
     }// end prime initialisation
   
     for i in 0..primes.len(){// sieve current segment
     
         let mut j = multiples[i as usize];
         
         let k = primes[i as usize]*2;
         
         loop {
           if j >= segment_size{
             break;
           }
           sieve[j as usize] = false;
           j+=k;
         }
         multiples[i as usize] = j -segment_size;
     }// end current sieve

     loop{
       if n > high{
         break;
       }
       if sieve[(n - low) as usize]{
         count+=1;
         
         if collect{
         values.push(n);
         }
       }
       n+=2;
     }
          low+=segment_size;
   }
   (values,count)
   }
