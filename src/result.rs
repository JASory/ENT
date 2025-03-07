/// Enum returned in checked functions returning any errors of evaluation 
#[non_exhaustive]
#[derive(Clone,Debug,PartialEq,Eq)]
pub enum NTResult<T: Sized + Clone>{
/// Evaluation is knowable and was calculated in the bounds. FFI return 0
  Eval(T), //  Flag 0
/// Solution Exists but does not fit in datatype. FFI return 1   
  Overflow,  //  Flag 1
/// Solution Does Not Exist. FFI return 2  
  DNE,       //  Flag 2
/// Solution Exists but is Infinitely large. FFI return 3   
  Infinite,  //  Flag 3
/// Solutions exist but there are infinite number of them (all elements of Z). FFI return 4
  InfiniteSet, // Flag 4
/// Computation exceeds some practical bound. FFI return 5
  CompExceeded,  //  Flag 5
/// Overflowed during computation result does not necessarily exceed the datatype. FFI return 6
  CompOverflow, // Flag 6
/// Function is not defined for the inputs. FFI return 7
  Undefined,  // Flag 7
}


impl<T: Sized + Clone + Default> NTResult<T>{

/// Returns the Evaluated number
 pub fn unwrap(&self) -> T{
    match self{
    NTResult::Eval(x) => x.clone(),
    NTResult::Overflow=> panic!("Result does not fit in datatype"),
    NTResult::DNE => panic!("Does Not Exist"),
    NTResult::Infinite => panic!("Infinitely large solution"),
    NTResult::InfiniteSet => panic!("Infinite number of solutions"),
    NTResult::CompExceeded => panic!("Computation exceeded hardcoded limit"),
    NTResult::CompOverflow => panic!("Overflowed during computation"),
    NTResult::Undefined => panic!("Undefined solution"),
    }
 }
 
 /// Returns the result or panics with the selected message
 pub fn expect(&self, signal: &str)-> T{
    match self{
     NTResult::Eval(x) => x.clone(),
     _=> panic!("{}",signal)
    }
 }
 /// Returns the result or the selected default
 pub fn unwrap_or(&self, res: T) -> T{
    match self{
     NTResult::Eval(x) => x.clone(),
     _=> res,
    }
 }
 /// Checks if the result is infinitely large
 pub fn is_infinite(&self) -> bool{
   match self{
   NTResult::Infinite => true,
     _=> false,
   }
 }
 /// Checks if the result is an infinite set
 pub fn is_infiniteset(&self) -> bool{
   match self{
   NTResult::InfiniteSet => true,
     _=> false,
   }
 }
 /// Checks if the solution to the function exists
 pub fn is_dne(&self) -> bool{
   match self{
   NTResult::DNE => true,
     _=> false,
   }
 }
 /// Checks if the result exceeds the datatype max
 pub fn is_overflow(&self) -> bool{
   match self{
   NTResult::Overflow => true,
     _=> false,
   }
 }
 /// Checks if there was an overflow during computation
 pub fn is_comp_overflow(&self) -> bool{
   match self{
   NTResult::CompOverflow => true,
     _=> false,
   }
 }
 /// Checks if the Computation exceeded 
 pub fn is_comp_exceeded(&self) -> bool{
   match self{
   NTResult::CompExceeded => true,
     _=> false,
   }
 }
  /// Checks if the result is undefined 
 pub fn is_undefined(&self) -> bool{
   match self{
   NTResult::Undefined => true,
     _=> false,
   }
 }
 /** Converts from Option to NTResult. None values are converted to the selected NTResult variant
 ```
  use number_theory::NumberTheory;
  use number_theory::NTResult; 
  
  // A None is returned here due to exceeding the i8::MAX
  let res = 255u8.checked_add(1);
  // The None option is converted to an NTResult::Overflow
   
  let convert = NTResult::from_option(res,NTResult::Overflow);
 
  assert_eq!(convert, NTResult::Overflow);
  
  // An Some result is converted to Eval 
  
    let res = 254u8.checked_add(1);
    let convert = NTResult::from_option(res,NTResult::Overflow);
    
      assert_eq!(convert, NTResult::Eval(255));
 ```
 */
 pub fn from_option( opt: Option<T>, ntvariant: Self) -> Self{
     match opt{
       Some(x) => NTResult::Eval(x),
       None => ntvariant,
     }
 }

/// Maps within the Enum
 pub fn map<U: Clone,F: FnOnce(T) -> U>(&self, func : F) -> NTResult<U>{
   match self {
     NTResult::Eval(x) => NTResult::Eval((func)(x.clone())),
     NTResult::Overflow => NTResult::Overflow,
     NTResult::DNE => NTResult::DNE,
     NTResult::Infinite => NTResult::Infinite,
     NTResult::InfiniteSet => NTResult::InfiniteSet,
     NTResult::CompExceeded => NTResult::CompExceeded,
     NTResult::CompOverflow => NTResult::CompOverflow,
     NTResult::Undefined => NTResult::Undefined,
   }
 }
 /** Return value and flag for FFI bindings
 ```
 use number_theory::NumberTheory;
   // Attempt to factor 0
   let res = 0.factor();
   // Fails and returns a NTResult that can be decomposed to an empty vec  
   // and 4 for C api binding
  
   
   let res = 1.factor(); 
   // Likewise an DNE NTResult gets converted to a vec and 2
 
   
   let res = 9.factor();
   // Finally a fully factorable integer gets a vector and the 0 flag
  
 ```  
 
 */ 
 pub fn ffi(&self)-> (T,u8){
    match self {
     NTResult::Eval(x) => (x.clone(),0u8),
     NTResult::Overflow => (T::default(),1u8),
     NTResult::DNE => (T::default(),2u8),
     NTResult::Infinite => (T::default(),3u8),
     NTResult::InfiniteSet => (T::default(),4u8),
     NTResult::CompExceeded => (T::default(),5u8),
     NTResult::CompOverflow => (T::default(),6u8),
     NTResult::Undefined => (T::default(),7u8),
   }
 }

}

impl<T: Clone + Default + Sized + std::fmt::Display > std::fmt::Display for NTResult<T>{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
       match self{
       
       NTResult::Eval(x) => write!(f,"{}",x),
       NTResult::Overflow => write!(f,"Overflow"),
       NTResult::DNE => write!(f,"DNE"),
       NTResult::Infinite => write!(f,"Infinite"),
       NTResult::InfiniteSet => write!(f,"Infinite solutions"),
       NTResult::CompExceeded => write!(f,"Exceeded computation bound"),
       NTResult::CompOverflow => write!(f,"Overflowed during computation"),
       NTResult::Undefined => write!(f,"Undefined"),
       }
    }
   } 
