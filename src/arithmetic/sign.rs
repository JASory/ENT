
/// Enum representing the sign
#[derive(PartialEq, Clone, Debug, Default)]
pub enum Sign {
   /// N >= 0 representation
    #[default]
    Positive,
   /// N <  0 representation 
    Negative,
}

impl Sign {
    pub(crate) fn neg(&self) -> Sign {
        match self {
            Sign::Positive => Sign::Negative,
            Sign::Negative => Sign::Positive,
        }
    }

    pub(crate) fn mul(&self, other: &Sign) -> Sign {
        match (self, other) {
            (&Sign::Positive, &Sign::Negative) => Sign::Negative,
            (&Sign::Negative, &Sign::Positive) => Sign::Negative,
            _ => Sign::Positive,
        }
    }

    pub(crate) fn _pow(&self, y: u64) -> Sign {
        if self == &Sign::Negative && y % 2 == 1 {
            Sign::Negative
        } else {
            Sign::Positive
        }
    }
    
 ///   To boolean for FFI
    pub fn ffi(&self) -> bool{
       match self{
         Sign::Negative => true,
         Sign::Positive => false,
       }
    }
  /// From boolean FFI  
    pub fn from_ffi(x: bool) -> Self{
    match x{
        true => Sign::Negative,
        false => Sign::Positive,
        }
    }
}
