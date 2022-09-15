It has been observed that the standard checks that are performed in most software are only sufficient for the most trivial errors, 
and in the field of mathematics detecting unknown errors involves much more statistical analysis than provided by quick edge case checks

As such number-theory (ENT) has two sets of checks trivial ones that are labeled "simple" and more analytic ones labeled "strong" that are 
vastly superior at detecting errors than the simple ones. Some of the strong checks may rigourously prove a function is correct, but more likely they give a statistical argument for it. It is the end goal of this library to be able to run checks that rigourously prove all functions to be correct. 

Note that these are not the extent of the checks performed by the author but rather the checks builtin for the user, the author has performed far stronger checks that
 are considered to take too long except for the most fastidious user (aka t > 1E+5 s)

{For the purposes of this library. "proven to be correct", means behaves exactly as described in the documentation}

List of functions fully proven to be correct via strong tests in this crate

is_prime      Proven to be deterministic under 2^64+2^42 and have an accuracy of atleast 2^-64 via Damgard et.al and heuristic check of generated pseudoprimes






Proven over intervals 










Broken functions 

logarithm(s)       Fails to account for  negative integers, which require a complex result, resolution can either be to panic or evaluate 

  
