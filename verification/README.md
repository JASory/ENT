
As ENT prioritizes correctness it follows a strict protocol. 

First we define a battery of tests that verify that the functions behave as described, the arguments may be statistical or proofs by exhaustion.
These tests are described in /verification/test and some are implemented in the libraries tests. However, due to the complexity and strength of the tests 
(handling gigabytes of data as in the case of checking that is_prime almost certainly does have a accuracy in excess of 2^-64) most are not. 

Iff (if and only if) the tests are all passed then we compute a pair of hashes for each source code file using Open-SSl's implementation of SHA-512-256 
and SHA-3-256 and write it into /verification/hash. This ensures that users can verify if the sourcecode does infact pass the strong tests by computing
the hash themselves and comparing it. Note that if a hash fails it does not necessarily mean that the sourcecode is broken, simply that it is not tested. 

As the tests have not been fully decided on, this protocol is not being strictly followed. Once ENT reaches 0.1.0 however this will be enforced for every
update. 
