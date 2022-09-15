/*

  List of computationally intensive tests that demonstrate the claims made by the author.


  Full computation of each individual test may take hours so it is imperative that it is run as cargo test --release and during downtime.

  These tests are primarily to ensure confidence in program correctness, especially in light of the errors presented by other number-theory libraries in Rust 3
  and elsewhere.

  The assumptions that are made to claim correctness are listed next to the test. If the assumptions are wrong then the proof of correctness is also wrong

  Assumption abbreviations :  P = mathematically proven, V = Verified by at least two parties including at least partially by the author,
   UV = Computed by the author, but unverified by external parties,
   SP = statistically probable, applies to analyzing external software that corroborates well with proven bounds

*/

use number_theory::Mpz;
use number_theory::NumberTheory;

/*

   Description                      Assumptions                                                   Approximate Time
_______________________________________________________________________________________________________________
 Exhaustive proof of            Correctness of Erastothenes sieve in computing pi(n) (P)              5E+3 s
 deterministic primality  	Correctness of strong fermat test (P)
 for integers under 		Correctness of Feitsma-Galway pseudoprime table (V)
 2^64 + 2^44                    Correctness of J.A Sory's reduction and extension of table (UV)
                                Correctness of Kim Walisch's implementation of Gourdon's pi(n) (SP)
                                Correctness of Damgard et .al in pseudoprime density (P)

_________________________________________________________________________________________________________________



*/

#[ignore]
#[test]
fn det_primality() {
    /* Determinism between 2^35;2^64 + 2^45  is easy as we only have to check the base-2 pseudoprimes
    and a table of such pseudoprimes has been provided by Jan Feitsma, William Galway.
    J.A Sory reduced the table to strong pseudoprimes coprime to the first 64 odd primes and extended it to 2^64+2^45.
    We use this variant below due to 4x faster evaluation.

    Contact j_a_sory at rust-cas dot org to obtain a copy as it is not currently publicly available.

    */

    let location = "/home/jasory/Proving/FGS2SPRP";
    let data = std::fs::read_to_string(location).expect(
        "Change your file path: Unable to find the base-2 pseudoprime table to prove correctness",
    );
    // no checks performed as all numbers are guaranteed to be valid if you are using Feitsma or Jasory's table
    let pseudos = data
        .split_whitespace()
        .map(|x| x.parse::<u128>().unwrap())
        .collect::<Vec<u128>>();

    for i in pseudos {
        assert_eq!(i.is_prime(), false)
    }

    let mut count = 0u32;

    for i in 0..u32::MAX {
        if i.is_prime() {
            count += 1;
        }
    }

    // Deterministic under 2^32, this must be computed exhaustively as it utilizes a single test
    assert_eq!(count, 203280221);

    // Finish checking up to 2^35
    for i in 4294967296u64..34359738368u64 {
        if i.is_prime() {
            count += 1;
        }
    }
    // Deterministic under 2^35
    assert_eq!(count, 1480206279)
}

/*
   Checks that approximately half (60%+) of the base-2 strong pseudoprimes between 2^64 and 2^67  are filtered. This is almost certainly an error rate much
   less than 2^-64.
   (Probabilistic claim as only about 1 base-2 pseudoprime exists per 2^40 composites in the interval. However fully computing pseudoprimes up to 2^67 is not possible,
   so this cannot be proven)

   Utilizes a pre-computed table of composites of the form (k+1)(rk+1) where r > 65535

   Fails at a rate of 5%, returning a 1, this is well-within the acceptable bounds for the advertised accuracy which is 4 per check (2^67/2^64 = 8*60% = 4)
*/
#[ignore]
#[test]
fn heuristic_test() {
    // List of approximately 60% of the 2-strong-pseudoprimes between 2^64;2^67 or 31280009 pseudoprimes
    let location = "/home/jasory/Proving/H2SPRP-64-67";
    let data = std::fs::read_to_string(location).expect(
        "Change your file path: Unable to find the base-2 pseudoprime table to prove correctness",
    );
    // no checks performed as all numbers are guaranteed to be valid
    let mut failures = 0u64;
    let h_pseudos = data
        .split_whitespace()
        .map(|x| x.parse::<u128>().unwrap())
        .collect::<Vec<u128>>();
    assert_eq!(h_pseudos.len(), 31280009);
    for i in h_pseudos {
        //assert_eq!(i.is_prime(), false)
        if i.is_prime() {
            failures += 1;
        }
    }
    assert_eq!(failures, 0) // may fail however as long as the result is less than 4 (never observed) the accuracy bound is valid
}

/*

 Tested                                               Permitted Error

 1E+3 iterations of Arnault's Carmichael number       30% (This is actually asymptotic to 1/4 however for smaller values it may exceed it so the permitted bound is 3/10)
 10 512+ bit Carmichael numbers                        0% // These Carmichael numbers are a subset of the form (6k+1)(12k+1)(18k+1)
 10 4000+ bit Carmichael numbers                       0%


  IS_Prime is possibly deterministic against Carmichael numbers of the form (6k+1)(12k+1)(18k+1) (This is hard to evaluate as none of )
*/

#[ignore]
#[test]
fn carmichael_test() {
    let mut failures = 0u64;

    let arnault = Mpz::u_new(vec![
        // Arnault's exceptionally strong Carmichael number that is strong pseudoprime to all n < 307, and approximately 1/3 of all bases
        14125993214864411435,
        17627421051211808665,
        8803062915072472141,
        16592776965248599063,
        2400842081300231610,
        15499637576960551660,
        13635884493322354475,
        11162303294968070407,
        16638949889857351534,
        4755447895431218293,
        3008559606663544016,
        6962752268970074012,
        15476370627519396605,
        7818823728409314089,
        7726208862163452462,
        2793947945622957183,
        9898516143788448943,
        2723958373070846137,
        15945921893246334102,
        9257474160255769327,
        138699416178,
    ]);

    for _ in 0..1000 {
        if arnault.is_prime() {
            failures += 1;
        }
    }

    assert!(failures < 300); // error rate is less than 30%.

    let six = Mpz::from_u64(6);
    let twelve = Mpz::from_u64(12);
    let eighteen = Mpz::from_u64(18);
    let one = Mpz::one();

    let rand_gen = |len: usize| -> (Mpz, Mpz, Mpz) {
        // Generates Carmichael numbers of the form (6k+1)(12k+1)(18k+1)
        loop {
            let rand = Mpz::rand(len, u64::rng);
            // rand.set_bit(0);

            let lhs = six.ref_product(&rand).ref_addition(&one);
            let mid = twelve.ref_product(&rand).ref_addition(&one);
            let rhs = eighteen.ref_product(&rand).ref_addition(&one);
            if lhs.is_prime() {
                if mid.is_prime() {
                    if rhs.is_prime() {
                        return (lhs, mid, rhs);
                    }
                }
            }
        }
    };

    for _ in 0..10 {
        // 512 bit test

        let (lhs, mid, rhs) = rand_gen(3);
        let carmichael = lhs.ref_product(&mid).ref_product(&rhs);
        assert_eq!(carmichael.is_prime(), false)
    }

    for _ in 0..10 {
        // 1E+4 Carmichael numbers of the same magnitude as Arnault's
        let (lhs, mid, rhs) = rand_gen(21);
        let carmichael = lhs.ref_product(&mid).ref_product(&rhs);
        assert_eq!(carmichael.is_prime(), false)
    }
}

/*

 Constructs 1E+5 probable semiprimes of the form (k+1)(rk+1) where r < 64. This is approximately half of the
  strong pseudoprimes. If this test passes then it is safe to judge the accuracy

  This test takes approximately 160 hrs to execute.

*/
#[ignore]
#[test]
fn probable_prime() {
    let rand_gen = |len: usize, k: u64| -> (Mpz, Mpz) {
        loop {
            let rng = Mpz::rand(len, u64::rng);
            let rk_one = Mpz::from_u64(k).ref_product(&rng).ref_addition(&Mpz::one());
            let k_one = rng.ref_addition(&Mpz::one());

            if k_one.is_prime() & rk_one.is_prime() {
                return (k_one, rk_one);
            }
        }
    };

    //let mut veccy = vec![];
    // generate and collect 10000 semiprimes of the form (k+1)(2k+1)
    let mut count = 0u64;
    for j in 2..64 {
        for _i in 0..100 {
            let (k_one, rk_one) = rand_gen(4, j);
            if k_one.ref_product(&rk_one).is_prime() {
                count += 1
            }
        }
    }
    assert_eq!(count, 0)
}

/* List of primes that are/were passed by prominent libraries as described in Prime & Prejudice: Primality Testing Albrecht et al .
Testing them here is largely irrelevant as the  purpose of the paper was  to demonstrate that deterministically selected bases fail
against selected primes and all deterministically selected bases in Number-Theory have either been proven to filter all composites or are complemented by random bases.

As implied this test is simply for user confidence rather than providing any particular insight.



*/

#[rustfmt::skip]
#[ignore]
#[test]

fn confidence_test(){

 let gmp_fear = Mpz::u_new(vec![ // Claimed to pass GMP 6.1.2 
        13347957477412617151, 16259308113041253372,
        8892814040333046162, 1221545748070377901, 
        13845990824863685200, 11427719796916674983, 
        6523989286207306293, 8228562809691301656, 
        15942565893161427335, 17798657415511144750, 
        9895556227828505149, 1067200624540728504, 6934903806235641717, 
        8868460805486687741, 8433825644186292052, 9355494833841550342
        ]); 
        
        for _ in 0..10000{// Deterministic due to a radical difference in base selection
         assert_eq!(gmp_fear.is_prime(), false)
        }
        
        
        // Claimed to deterministically pass Apple, Libtom, CommonCrypto & WolfSSL, see original paper for software versions (this may no longer be true);
        // Note that for composites of the form p*(x(p-1)+1)(y(p-1)+1), is_prime has at worst around 1/4 
  let security_flaw =     Mpz::u_new(vec![
  9331929281010904555, 5494805439598348446, 8274855319168077407, 
8989115052798713304, 7922897861797270742, 16684056220225637685,
11573383658562860865, 10538110707526527084, 11841525957671433530,
6629998598229579978, 9754395927480723018, 1082083247965744052,
13178505664201790977, 445352704338639786, 1262458056955046950, 
10952769848018447964, 9675180056116935648, 18095375587926220513, 
14906862786560135370, 7349764202741105504, 14310555118399348948, 
16532416774430072831, 5249926442601408071, 10781540968865591231, 
15478470103784080673, 12574244713382308730, 7937862675332072157, 
8979440646405107343, 7057598508513752002, 330859835934063533,
15284532345271639087, 2683089657284425004, 12388650678928717386, 
4601150668493820598, 17591554493017240341, 8874420574259649282,
5916665853078007903, 1744435599167408248, 1375564577028398640,
11421316388837429873, 16721851356632668685, 3143746952704220612, 
9964885490987307667, 17847312369004334358, 2102167341425916716, 
17402799351753590688, 5031694375487631134, 9769474730685880473, 
9334559868889485670, 5705824887076732981, 10296763285808252259, 
15052706519152474971, 10441944543730280361, 3742633571041123035, 
7671374831153001193, 5651690108603415721, 7360676820152831689, 
3300584806060773119, 13951650231935285153, 14432958795279723086, 
5197204459193853662, 16439386015664755780, 15489184032414823641, 
5445353960443525377, 17365100034401419080, 2507657222167183359, 
7305097783839097981, 3577408235263180558, 17574691796078419858, 
7220852877782605389, 5724975390124524347, 7983361366552541545, 
7183619713474561234, 14313726161713768561, 17635332671748993665, 
10189608807876133269, 15326339027698575730, 16262806750666400113, 
6718975974224483036, 3265280417026005297, 8262935061937412279, 
11470300983668961140, 18234435204737306426, 13868574488010990731, 
15544778790353478078, 1654287122911138553, 12115954827342798448, 
18315805652863710573, 9072519746302664175, 9530505032704430404, 
10950413179334519245, 11500900080207907576, 5289500343761324225, 
5730607270653974271, 12934085552529239836, 1627559471621096679, 
330020620966697801, 15559778274953646736, 11080068881085564969, 
10399807091365952939, 4843032489303663269, 16516455478232626003, 
10601252522603802158, 767100916292225657, 2497006575850938970, 
5676095330871291613, 11532840839212411388, 14210693347085223154, 
13738501162309476745, 110019121374055]
  );  
        let mut failures = 0u8;
     for _ in 0..100{
         if security_flaw.is_prime(){
           failures+=1
         }
        }
  assert!(failures < 30)
}
