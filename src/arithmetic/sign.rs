#[derive(PartialEq,Clone, Debug)]
pub enum Sign{
                Positive,
                Negative,
}

impl Sign{

pub(crate) fn neg(&self) -> Sign {
              match self {
              Sign::Positive =>  Sign::Negative,
              Sign::Negative =>  Sign::Positive, 
              }
           }

pub(crate) fn mul(&self, other: &Sign)-> Sign {
             match (self, other){
              (&Sign::Positive, &Sign::Negative)=>  Sign::Negative,
              (&Sign::Negative, &Sign::Positive)=>  Sign::Negative,
                                        _=>  Sign::Positive,
            }
 }
 
pub(crate) fn pow(&self, y: u64) -> Sign {
               if self == &Sign::Negative && y%2 ==1{
                  return Sign::Negative
               }
               else{
                 Sign::Positive
               }
              }  
  }
