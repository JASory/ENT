/*

      Prime Data used for primality testing




  "DET_MAX" is the upperbound for deterministic tests, all primes below this number receive a maximum of 2 sprp tests,
   a considerable speed up over the minimum of 12 tests that have been previously proven, note that
   this is approximately  145 trillion higher than the bound of 2^64, (and over 3 trillion more primes) provided by other tests and
   continously increases due to research by J.A Sory
*/
pub(crate) const DET_MAX: u128 = 0x10000840000000000; // current bound 2^64 + 2^47 + 2^42

// List of Mersenne prime exponents, shortcuts computation as all Mersenne's not listed below the bound of 57885161 have been proven composite
pub(crate) const MERSENNE_LIST: [u32; 42] = [
    89, 107, 127, 521, 607, 1279, 2203, 2281, 3217, 4253, 4423, 9689, 9941, 11213, 19937, 21701,
    23209, 44497, 86243, 110503, 132049, 216091, 756839, 859433, 1257787, 1398269, 2976221,
    3021377, 6972593, 13466917, 20996011, 24036583, 25964951, 30402457, 32582657, 37156667,
    42643801, 43112609, 57885161, 74207281, 77232917, 82589933,
]; 

// Multiplicative inverse of odd integers in Z[2^8]
#[rustfmt::skip]
pub(crate) const INV_8 : [u8; 128] = [

    0x01, 0xAB, 0xCD, 0xB7, 0x39, 0xA3, 0xC5, 0xEF, 0xF1, 0x1B, 0x3D, 0xA7, 0x29, 0x13, 0x35, 0xDF,
    0xE1, 0x8B, 0xAD, 0x97, 0x19, 0x83, 0xA5, 0xCF, 0xD1, 0xFB, 0x1D, 0x87, 0x09, 0xF3, 0x15, 0xBF,
    0xC1, 0x6B, 0x8D, 0x77, 0xF9, 0x63, 0x85, 0xAF, 0xB1, 0xDB, 0xFD, 0x67, 0xE9, 0xD3, 0xF5, 0x9F,
    0xA1, 0x4B, 0x6D, 0x57, 0xD9, 0x43, 0x65, 0x8F, 0x91, 0xBB, 0xDD, 0x47, 0xC9, 0xB3, 0xD5, 0x7F,
    0x81, 0x2B, 0x4D, 0x37, 0xB9, 0x23, 0x45, 0x6F, 0x71, 0x9B, 0xBD, 0x27, 0xA9, 0x93, 0xB5, 0x5F,
    0x61, 0x0B, 0x2D, 0x17, 0x99, 0x03, 0x25, 0x4F, 0x51, 0x7B, 0x9D, 0x07, 0x89, 0x73, 0x95, 0x3F,
    0x41, 0xEB, 0x0D, 0xF7, 0x79, 0xE3, 0x05, 0x2F, 0x31, 0x5B, 0x7D, 0xE7, 0x69, 0x53, 0x75, 0x1F,
    0x21, 0xCB, 0xED, 0xD7, 0x59, 0xC3, 0xE5, 0x0F, 0x11, 0x3B, 0x5D, 0xC7, 0x49, 0x33, 0x55, 0xFF,

];

// Multiplicative inverse of the first 256 odd primes (3..1621 inclusive), in Z[2^64]
#[rustfmt::skip]
pub(crate) static PRIME_INV_64: [u64; 256] = [
    
           0xaaaaaaaaaaaaaaab ,   0xcccccccccccccccd ,   0x6db6db6db6db6db7 ,   0x2e8ba2e8ba2e8ba3 ,
	   0x4ec4ec4ec4ec4ec5 ,   0xf0f0f0f0f0f0f0f1 ,   0x86bca1af286bca1b ,   0xd37a6f4de9bd37a7 ,
	   0x34f72c234f72c235 ,   0xef7bdef7bdef7bdf ,   0x14c1bacf914c1bad ,   0x8f9c18f9c18f9c19 ,
	   0x82fa0be82fa0be83 ,   0x51b3bea3677d46cf ,   0x21cfb2b78c13521d ,   0xcbeea4e1a08ad8f3 ,
	   0x4fbcda3ac10c9715 ,   0xf0b7672a07a44c6b ,   0x193d4bb7e327a977 ,   0x7e3f1f8fc7e3f1f9 ,
	   0x9b8b577e613716af ,   0xa3784a062b2e43db ,   0xf47e8fd1fa3f47e9 ,   0xa3a0fd5c5f02a3a1 ,
	   0x3a4c0a237c32b16d ,   0xdab7ec1dd3431b57 ,   0x77a04c8f8d28ac43 ,   0xa6c0964fda6c0965 ,
	   0x90fdbc090fdbc091 ,   0x7efdfbf7efdfbf7f ,    0x3e88cb3c9484e2b ,   0xe21a291c077975b9 ,
	   0x3aef6ca970586723 ,   0xdf5b0f768ce2cabd ,   0x6fe4dfc9bf937f27 ,   0x5b4fe5e92c0685b5 ,
	   0x1f693a1c451ab30b ,   0x8d07aa27db35a717 ,   0x882383b30d516325 ,   0xed6866f8d962ae7b ,
	   0x3454dca410f8ed9d ,   0x1d7ca632ee936f3f ,   0x70bf015390948f41 ,   0xc96bdb9d3d137e0d ,
	   0x2697cc8aef46c0f7 ,   0xc0e8f2a76e68575b ,   0x687763dfdb43bb1f ,   0x1b10ea929ba144cb ,
	   0x1d10c4c0478bbced ,   0x63fb9aeb1fdcd759 ,   0x64afaa4f437b2e0f ,   0xf010fef010fef011 ,
	   0x28cbfbeb9a020a33 ,   0xff00ff00ff00ff01 ,   0xd624fd1470e99cb7 ,   0x8fb3ddbd6205b5c5 ,
	   0xd57da36ca27acdef ,   0xee70c03b25e4463d ,   0xc5b1a6b80749cb29 ,   0x47768073c9b97113 ,
	   0x2591e94884ce32ad ,   0xf02806abc74be1fb ,   0x7ec3e8f3a7198487 ,   0x58550f8a39409d09 ,
	   0xec9e48ae6f71de15 ,   0x2ff3a018bfce8063 ,   0x7f9ec3fcf61fe7b1 ,   0x89f5abe570e046d3 ,
	   0xda971b23f1545af5 ,   0x79d5f00b9a7862a1 ,   0x4dba1df32a128a57 ,   0x87530217b7747d8f ,
	   0x30baae53bb5e06dd ,   0xee70206c12e9b5b3 ,   0xcdde9462ec9dbe7f ,   0xafb64b05ec41cf4d ,
	    0x2944ff5aec02945 ,   0x2cb033128382df71 ,   0x1ccacc0c84b1c2a9 ,   0x19a93db575eb3a0b ,
	   0xcebeef94fa86fe2d ,   0x6faa77fb3f8df54f ,   0x68a58af00975a751 ,   0xd56e36d0c3efac07 ,
	   0xd8b44c47a8299b73 ,    0x2d9ccaf9ba70e41 ,    0x985e1c023d9e879 ,   0x2a343316c494d305 ,
	   0x70cb7916ab67652f ,   0xd398f132fb10fe5b ,   0x6f2a38a6bf54fa1f ,   0x211df689b98f81d7 ,
	    0xe994983e90f1ec3 ,   0xad671e44bed87f3b ,   0xf9623a0516e70fc7 ,   0x4b7129be9dece355 ,
	   0x190f3b7473f62c39 ,   0x63dacc9aad46f9a3 ,   0xc1108fda24e8d035 ,   0xb77578472319bd8b ,
	   0x473d20a1c7ed9da5 ,   0xfbe85af0fea2c8fb ,   0x58a1f7e6ce0f4c09 ,   0x1a00e58c544986f3 ,
	   0x7194a17f55a10dc1 ,   0x7084944785e33763 ,   0xba10679bd84886b1 ,   0xebe9c6bb31260967 ,
	   0x97a3fe4bd1ff25e9 ,   0x6c6388395b84d99f ,   0x8c51da6a1335df6d ,   0x46f3234475d5add9 ,
	   0x905605ca3c619a43 ,   0xcee8dff304767747 ,   0xff99c27f00663d81 ,   0xacca407f671ddc2b ,
	   0xe71298bac1e12337 ,   0xfa1e94309cd09045 ,   0xbebccb8e91496b9b ,   0x312fa30cc7d7b8bd ,
	   0x6160ff9e9f006161 ,   0x6b03673b5e28152d ,   0xfe802ffa00bfe803 ,   0xe66fe25c9e907c7b ,
	   0x3f8b236c76528895 ,   0xf6f923bf01ce2c0d ,   0x6c3d3d98bed7c42f ,   0x30981efcd4b010e7 ,
	   0x6f691fc81ebbe575 ,   0xb10480ddb47b52cb ,   0x74cd59ed64f3f0d7 ,    0x105cb81316d6c0f ,
	   0x9be64c6d91c1195d ,   0x71b3f945a27b1f49 ,   0x77d80d50e508fd01 ,   0xa5eb778e133551cd ,
	   0x18657d3c2d8a3f1b ,   0x2e40e220c34ad735 ,   0xa76593c70a714919 ,   0x1eef452124eea383 ,
	   0x38206dc242ba771d ,   0x4cd4c35807772287 ,   0x83de917d5e69ddf3 ,   0x882ef0403b4a6c15 ,
	   0xf8fb6c51c606b677 ,   0xb4abaac446d3e1fd ,   0xa9f83bbe484a14e9 ,    0xbebbc0d1ce874d3 ,
	   0xbd418eaf0473189f ,   0x44e3af6f372b7e65 ,   0xc87fdace4f9e5d91 ,   0xec93479c446bd9bb ,
	   0xdac4d592e777c647 ,   0xa63ea8c8f61f0c23 ,   0xe476062ea5cbbb6f ,   0xdf68761c69daac27 ,
	   0xb813d737637aa061 ,   0xa3a77aac1fb15099 ,   0x17f0c3e0712c5825 ,   0xfd912a70ff30637b ,
	   0xfbb3b5dc01131289 ,   0x856d560a0f5acdf7 ,   0x96472f314d3f89e3 ,   0xa76f5c7ed2253531 ,
	   0x816eae7c7bf69fe7 ,   0xb6a2bea4cfb1781f ,   0xa3900c53318e81ed ,   0x60aa7f5d9f148d11 ,
	   0x6be8c0102c7a505d ,   0x8ff3f0ed28728f33 ,   0x680e0a87e5ec7155 ,   0xbbf70fa49fe829b7 ,
	   0xd69d1e7b6a50ca39 ,   0x1a1e0f46b6d26aef ,   0x7429f9a7a8251829 ,   0xd9c2219d1b863613 ,
	   0x91406c1820d077ad ,   0x521f4ec02e3d2b97 ,   0xbb8283b63dc8eba5 ,   0x431eda153229ebbf ,
	   0xaf0bf78d7e01686b ,   0xa9ced0742c086e8d ,   0xc26458ad9f632df9 ,   0xbbff1255dff892af ,
	   0xcbd49a333f04d8fd ,   0xec84ed6f9cfdeff5 ,   0x97980cc40bda9d4b ,   0x777f34d524f5cbd9 ,
	   0x2797051d94cbbb7f ,   0xea769051b4f43b81 ,   0xce7910f3034d4323 ,   0x92791d1374f5b99b ,
	   0x89a5645cc68ea1b5 ,   0x5f8aacf796c0cf0b ,   0xf2e90a15e33edf99 ,   0x8e99e5feb897c451 ,
	   0xaca2eda38fb91695 ,   0x5d9b737be5ea8b41 ,   0x4aefe1db93fd7cf7 ,   0xa0994ef20b3f8805 ,
	   0x103890bda912822f ,   0xb441659d13a9147d ,   0x1e2134440c4c3f21 ,   0x263a27727a6883c3 ,
	   0x78e221472ab33855 ,   0x95eac88e82e6faff ,   0xf66c258317be8dab ,    0x9ee202c7cb91939 ,
	   0x8d2fca1042a09ea3 ,   0x82779c856d8b8bf1 ,   0x3879361cba8a223d ,   0xf23f43639c3182a7 ,
	   0xa03868fc474bcd13 ,   0x651e78b8c5311a97 ,   0x8ffce639c00c6719 ,   0xf7b460754b0b61cf ,
	   0x7b03f3359b8e63b1 ,   0xa55c5326041eb667 ,   0x647f88ab896a76f5 ,   0x8fd971434a55a46d ,
	   0x9fbf969958046447 ,   0x9986feba69be3a81 ,   0xa668b3e6d053796f ,   0x97694e6589f4e09b ,
	   0x37890c00b7721dbd ,   0x5ac094a235f37ea9 ,   0x31cff775f2d5d65f ,   0xddad8e6b36505217 ,
	   0x5a27df897062cd03 ,   0xe2396fe0fdb5a625 ,   0xb352a4957e82317b ,   0xd8ab3f2c60c2ea3f ,
	   0x6893f702f0452479 ,   0x9686fdc182acf7e3 ,   0x6854037173dce12f ,   0x7f0ded1685c27331 ,
	   0xeeda72e1fe490b7d ,   0x9e7bfc959a8e6e53 ,   0x49b314d6d4753dd7 ,   0x2e8f8c5ac4aa1b3b ,
	   0xb8ef723481163d33 ,   0x6a2ec96a594287b7 ,   0xdba41c6d13aab8c5 ,   0xc2adbe648dc3aaf1 ,
	   0x87a2bade565f91a7 ,   0x4d6fe8798c01f5df ,   0x3791310c8c23d98b ,   0xf80e446b01228883 ,
	   0x9aed1436fbf500cf ,   0x7839b54cc8b24115 ,   0xc128c646ad0309c1 ,   0x14de631624a3c377 ,
	   0x3f7b9fe68b0ecbf9 ,   0x284ffd75ec00a285 ,   0x37803cb80dea2ddb ,   0x86b63f7c9ac4c6fd ,
	
];

#[rustfmt::skip]
pub(crate) static PRIME_INV_128: [u128; 128] = [
    
          0xaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaab ,  0xcccccccccccccccccccccccccccccccd ,  
          0xb6db6db6db6db6db6db6db6db6db6db7 ,  0xa2e8ba2e8ba2e8ba2e8ba2e8ba2e8ba3 ,
	  0xc4ec4ec4ec4ec4ec4ec4ec4ec4ec4ec5 ,  0xf0f0f0f0f0f0f0f0f0f0f0f0f0f0f0f1 ,  
	  0xbca1af286bca1af286bca1af286bca1b ,  0x4de9bd37a6f4de9bd37a6f4de9bd37a7 ,
	  0xc234f72c234f72c234f72c234f72c235 ,  0xdef7bdef7bdef7bdef7bdef7bdef7bdf ,  
	  0xc1bacf914c1bacf914c1bacf914c1bad ,  0x18f9c18f9c18f9c18f9c18f9c18f9c19 ,
	  0xbe82fa0be82fa0be82fa0be82fa0be83 ,  0x3677d46cefa8d9df51b3bea3677d46cf ,  
	  0x13521cfb2b78c13521cfb2b78c13521d ,  0x8f2fba9386822b63cbeea4e1a08ad8f3 ,
	  0x14fbcda3ac10c9714fbcda3ac10c9715 ,  0xc2dd9ca81e9131abf0b7672a07a44c6b ,  
	  0x4f52edf8c9ea5dbf193d4bb7e327a977 ,  0x3f1f8fc7e3f1f8fc7e3f1f8fc7e3f1f9 ,
	  0xd5df984dc5abbf309b8b577e613716af ,  0x2818acb90f6bf3a9a3784a062b2e43db ,  
	  0xd1fa3f47e8fd1fa3f47e8fd1fa3f47e9 ,  0x5f02a3a0fd5c5f02a3a0fd5c5f02a3a1 ,
	  0xc32b16cfd7720f353a4c0a237c32b16d ,  0xd0c6d5bf60ee9a18dab7ec1dd3431b57 ,  
	  0xa2b10bf66e0e5aea77a04c8f8d28ac43 ,  0xc0964fda6c0964fda6c0964fda6c0965 ,
	  0xc090fdbc090fdbc090fdbc090fdbc091 ,  0xbf7efdfbf7efdfbf7efdfbf7efdfbf7f ,  
	  0xf82ee6986d6f63aa03e88cb3c9484e2b ,  0x21a291c077975b8fe21a291c077975b9 ,
	  0xa2126ad1f4f31ba03aef6ca970586723 ,  0x93c225cc74d50c06df5b0f768ce2cabd ,  
	  0x26fe4dfc9bf937f26fe4dfc9bf937f27 ,   0x685b4fe5e92c0685b4fe5e92c0685b5 ,
	  0x8bc775ca99ea03241f693a1c451ab30b ,  0x513ed9ad38b7f3bc8d07aa27db35a717 ,  
	  0x133caba736c05eb4882383b30d516325 ,   0xe4d3aa30a02dc3eed6866f8d962ae7b ,
	  0x6fbc1c498c05a84f3454dca410f8ed9d ,  0x7749b79f7f5470961d7ca632ee936f3f ,  
	  0x90948f40feac6f6b70bf015390948f41 ,   0xbb207cc0532ae21c96bdb9d3d137e0d ,
	  0x7a3607b7f5b5630e2697cc8aef46c0f7 ,  0x2f514a026d31be7bc0e8f2a76e68575b ,  
	  0xdd8f7f6d0eec7bfb687763dfdb43bb1f ,  0x766a024168e18cf81b10ea929ba144cb ,
	   0xc4c0478bbcecfee1d10c4c0478bbced ,  0x758fee6bac7f735d63fb9aeb1fdcd759 ,   
	   0x77f76e538c5167e64afaa4f437b2e0f ,  0x10fef010fef010fef010fef010fef011 ,
	  0xa020a32fefae680828cbfbeb9a020a33 ,  0xff00ff00ff00ff00ff00ff00ff00ff01 ,  
	  0xf836826ef73d52bcd624fd1470e99cb7 ,  0x3ce8354b2ea1c8cd8fb3ddbd6205b5c5 ,
	  0x8715ba188f963302d57da36ca27acdef ,  0xb25e4463cff13686ee70c03b25e4463d ,  
	  0x6c69ae01d272ca3fc5b1a6b80749cb29 ,  0xf26e5c44bfc61b2347768073c9b97113 ,
	  0xb07dd0d1b15d7cf12591e94884ce32ad ,  0xd2f87ebfcaa1c5a0f02806abc74be1fb ,  
	  0xbe25dd6d7aa646ca7ec3e8f3a7198487 ,  0xbc1d71afd8bdc03458550f8a39409d09 ,
	  0x2ed6d05a72acd1f7ec9e48ae6f71de15 ,  0x62ff3a018bfce8062ff3a018bfce8063 ,  
	  0x3fcf61fe7b0ff3d87f9ec3fcf61fe7b1 ,  0x398b6f668c2c43df89f5abe570e046d3 ,
	  0x8c1a682913ce1eceda971b23f1545af5 ,   0xb9a7862a0ff465879d5f00b9a7862a1 ,  
	  0xe7c13f77161b18f54dba1df32a128a57 ,  0x73186a06f9b8d9a287530217b7747d8f ,
	  0x7c39a6c708ec18b530baae53bb5e06dd ,  0x37634af9ebbc742dee70206c12e9b5b3 ,  
	  0x503578fb5236cf34cdde9462ec9dbe7f ,  0xbcdfc0d2975ccab1afb64b05ec41cf4d ,
	  0xf5aec02944ff5aec02944ff5aec02945 ,  0xc7d208f00a36e71a2cb033128382df71 ,  
	  0xd38f55c0280f05a21ccacc0c84b1c2a9 ,  0xca3be03aa7687a3219a93db575eb3a0b ,
	  0x6a69ce2344b66c3ccebeef94fa86fe2d ,  0xfecfe37d53bfd9fc6faa77fb3f8df54f ,  
	  0xa58af00975a750ff68a58af00975a751 ,  0xdc6da187df580dfed56e36d0c3efac07 ,
	  0x8fe44308ab0d4a8bd8b44c47a8299b73 ,  0xf1bf0091f5bcb8bb02d9ccaf9ba70e41 ,  
	  0x5e1c023d9e878ff70985e1c023d9e879 ,  0x7880d53da3d15a842a343316c494d305 ,
	  0x1ddb81ef699b5e8c70cb7916ab67652f ,  0xf364512170607acad398f132fb10fe5b ,  
	  0xadb1f8848af4c6d06f2a38a6bf54fa1f ,  0xd9a0541b55af0c17211df689b98f81d7 ,
	  0x673bf5928258a2ac0e994983e90f1ec3 ,   0xdda093c0628041aad671e44bed87f3b ,  
	  0xa9fcf24229bbcd1af9623a0516e70fc7 ,  0xcbb18a4f7732cc324b7129be9dece355 ,
	   0x1f727cce5f530a5190f3b7473f62c39 ,  0x6da4f4bdeb71121c63dacc9aad46f9a3 ,  
	  0x4d9abc552cf42b88c1108fda24e8d035 ,  0x141fd3124095c328b77578472319bd8b ,
	  0xddfd3e0bf3218d19473d20a1c7ed9da5 ,  0xdb2b3278f3b910d2fbe85af0fea2c8fb ,  
	  0xcb5c3b636e3a7d1358a1f7e6ce0f4c09 ,  0x1bcbfe34e7576cf21a00e58c544986f3 ,
	  0x6b5e80aa5ef23f007194a17f55a10dc1 ,  0x9a628feb11022e3a7084944785e33763 ,  
	  0xbe61909edde53c01ba10679bd84886b1 ,  0x4feb7c5e05fbb9e8ebe9c6bb31260967 ,
	  0x1ff25e8ff92f47fc97a3fe4bd1ff25e9 ,  0x30143e6b1fa187616c6388395b84d99f ,  
	  0xd49154c6c94ac0f08c51da6a1335df6d ,  0x9b9771454a44e00d46f3234475d5add9 ,
	  0x3aba1b4baef0b2a9905605ca3c619a43 ,  0xcc11d9dd1bfe608ecee8dff304767747 ,  
	  0xff99c27f00663d80ff99c27f00663d81 ,  0x111ea8032f60bf1aacca407f671ddc2b ,
	  0xdd9395f5b667aa88e71298bac1e12337 ,  0xa7caaed93038740afa1e94309cd09045 ,  
	  0x2be5958f582e9db7bebccb8e91496b9b ,  0x995e1ca8dbfb5a3d312fa30cc7d7b8bd ,
	  0x9f006160ff9e9f006160ff9e9f006161 ,  0xb33ce15ee9b097416b03673b5e28152d ,  
	  0xfa00bfe802ffa00bfe802ffa00bfe803 ,  0x1c2802f6bcf18d26e66fe25c9e907c7b ,
	  0xcf6dec4793e72aba3f8b236c76528895 ,  0x1e547da72d224d44f6f923bf01ce2c0d ,  
	  0x7746da9d5fc708306c3d3d98bed7c42f ,  0xcdff4bb55916e37a30981efcd4b010e7 ,
	  
	 
];

// List of 2048 primes, exclusively used for trial division of huge integers, machine-sized words use precomputed inverses, Possibly remove?
pub(crate) const PRIMELIST :  [u16;2048] =  [
   2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
   73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
   179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281,
   283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409,
   419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541,
   547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659,
   661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809,
   811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941,
   947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013, 1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069,
   1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, 1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223,
   1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, 1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373,
   1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, 1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511,
   1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, 1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657,
   1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, 1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811,
   1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889, 1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987,
   1993, 1997, 1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053, 2063, 2069, 2081, 2083, 2087, 2089, 2099, 2111, 2113, 2129,
   2131, 2137, 2141, 2143, 2153, 2161, 2179, 2203, 2207, 2213, 2221, 2237, 2239, 2243, 2251, 2267, 2269, 2273, 2281, 2287,
   2293, 2297, 2309, 2311, 2333, 2339, 2341, 2347, 2351, 2357, 2371, 2377, 2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423,
   2437, 2441, 2447, 2459, 2467, 2473, 2477, 2503, 2521, 2531, 2539, 2543, 2549, 2551, 2557, 2579, 2591, 2593, 2609, 2617,
   2621, 2633, 2647, 2657, 2659, 2663, 2671, 2677, 2683, 2687, 2689, 2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731, 2741,
   2749, 2753, 2767, 2777, 2789, 2791, 2797, 2801, 2803, 2819, 2833, 2837, 2843, 2851, 2857, 2861, 2879, 2887, 2897, 2903,
   2909, 2917, 2927, 2939, 2953, 2957, 2963, 2969, 2971, 2999, 3001, 3011, 3019, 3023, 3037, 3041, 3049, 3061, 3067, 3079,
   3083, 3089, 3109, 3119, 3121, 3137, 3163, 3167, 3169, 3181, 3187, 3191, 3203, 3209, 3217, 3221, 3229, 3251, 3253, 3257,
   3259, 3271, 3299, 3301, 3307, 3313, 3319, 3323, 3329, 3331, 3343, 3347, 3359, 3361, 3371, 3373, 3389, 3391, 3407, 3413,
   3433, 3449, 3457, 3461, 3463, 3467, 3469, 3491, 3499, 3511, 3517, 3527, 3529, 3533, 3539, 3541, 3547, 3557, 3559, 3571,
   3581, 3583, 3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643, 3659, 3671, 3673, 3677, 3691, 3697, 3701, 3709, 3719, 3727,
   3733, 3739, 3761, 3767, 3769, 3779, 3793, 3797, 3803, 3821, 3823, 3833, 3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907,
   3911, 3917, 3919, 3923, 3929, 3931, 3943, 3947, 3967, 3989, 4001, 4003, 4007, 4013, 4019, 4021, 4027, 4049, 4051, 4057,
   4073, 4079, 4091, 4093, 4099, 4111, 4127, 4129, 4133, 4139, 4153, 4157, 4159, 4177, 4201, 4211, 4217, 4219, 4229, 4231,
   4241, 4243, 4253, 4259, 4261, 4271, 4273, 4283, 4289, 4297, 4327, 4337, 4339, 4349, 4357, 4363, 4373, 4391, 4397, 4409,
   4421, 4423, 4441, 4447, 4451, 4457, 4463, 4481, 4483, 4493, 4507, 4513, 4517, 4519, 4523, 4547, 4549, 4561, 4567, 4583,
   4591, 4597, 4603, 4621, 4637, 4639, 4643, 4649, 4651, 4657, 4663, 4673, 4679, 4691, 4703, 4721, 4723, 4729, 4733, 4751,
   4759, 4783, 4787, 4789, 4793, 4799, 4801, 4813, 4817, 4831, 4861, 4871, 4877, 4889, 4903, 4909, 4919, 4931, 4933, 4937,
   4943, 4951, 4957, 4967, 4969, 4973, 4987, 4993, 4999, 5003, 5009, 5011, 5021, 5023, 5039, 5051, 5059, 5077, 5081, 5087,
   5099, 5101, 5107, 5113, 5119, 5147, 5153, 5167, 5171, 5179, 5189, 5197, 5209, 5227, 5231, 5233, 5237, 5261, 5273, 5279,
   5281, 5297, 5303, 5309, 5323, 5333, 5347, 5351, 5381, 5387, 5393, 5399, 5407, 5413, 5417, 5419, 5431, 5437, 5441, 5443,
   5449, 5471, 5477, 5479, 5483, 5501, 5503, 5507, 5519, 5521, 5527, 5531, 5557, 5563, 5569, 5573, 5581, 5591, 5623, 5639,
   5641, 5647, 5651, 5653, 5657, 5659, 5669, 5683, 5689, 5693, 5701, 5711, 5717, 5737, 5741, 5743, 5749, 5779, 5783, 5791,
   5801, 5807, 5813, 5821, 5827, 5839, 5843, 5849, 5851, 5857, 5861, 5867, 5869, 5879, 5881, 5897, 5903, 5923, 5927, 5939,
   5953, 5981, 5987, 6007, 6011, 6029, 6037, 6043, 6047, 6053, 6067, 6073, 6079, 6089, 6091, 6101, 6113, 6121, 6131, 6133,
   6143, 6151, 6163, 6173, 6197, 6199, 6203, 6211, 6217, 6221, 6229, 6247, 6257, 6263, 6269, 6271, 6277, 6287, 6299, 6301,
   6311, 6317, 6323, 6329, 6337, 6343, 6353, 6359, 6361, 6367, 6373, 6379, 6389, 6397, 6421, 6427, 6449, 6451, 6469, 6473,
   6481, 6491, 6521, 6529, 6547, 6551, 6553, 6563, 6569, 6571, 6577, 6581, 6599, 6607, 6619, 6637, 6653, 6659, 6661, 6673,
   6679, 6689, 6691, 6701, 6703, 6709, 6719, 6733, 6737, 6761, 6763, 6779, 6781, 6791, 6793, 6803, 6823, 6827, 6829, 6833,
   6841, 6857, 6863, 6869, 6871, 6883, 6899, 6907, 6911, 6917, 6947, 6949, 6959, 6961, 6967, 6971, 6977, 6983, 6991, 6997,
   7001, 7013, 7019, 7027, 7039, 7043, 7057, 7069, 7079, 7103, 7109, 7121, 7127, 7129, 7151, 7159, 7177, 7187, 7193, 7207,
   7211, 7213, 7219, 7229, 7237, 7243, 7247, 7253, 7283, 7297, 7307, 7309, 7321, 7331, 7333, 7349, 7351, 7369, 7393, 7411,
   7417, 7433, 7451, 7457, 7459, 7477, 7481, 7487, 7489, 7499, 7507, 7517, 7523, 7529, 7537, 7541, 7547, 7549, 7559, 7561,
   7573, 7577, 7583, 7589, 7591, 7603, 7607, 7621, 7639, 7643, 7649, 7669, 7673, 7681, 7687, 7691, 7699, 7703, 7717, 7723,
   7727, 7741, 7753, 7757, 7759, 7789, 7793, 7817, 7823, 7829, 7841, 7853, 7867, 7873, 7877, 7879, 7883, 7901, 7907, 7919,
   7927, 7933, 7937, 7949, 7951, 7963, 7993, 8009, 8011, 8017, 8039, 8053, 8059, 8069, 8081, 8087, 8089, 8093, 8101, 8111,
   8117, 8123, 8147, 8161, 8167, 8171, 8179, 8191, 8209, 8219, 8221, 8231, 8233, 8237, 8243, 8263, 8269, 8273, 8287, 8291,
   8293, 8297, 8311, 8317, 8329, 8353, 8363, 8369, 8377, 8387, 8389, 8419, 8423, 8429, 8431, 8443, 8447, 8461, 8467, 8501,
   8513, 8521, 8527, 8537, 8539, 8543, 8563, 8573, 8581, 8597, 8599, 8609, 8623, 8627, 8629, 8641, 8647, 8663, 8669, 8677,
   8681, 8689, 8693, 8699, 8707, 8713, 8719, 8731, 8737, 8741, 8747, 8753, 8761, 8779, 8783, 8803, 8807, 8819, 8821, 8831,
   8837, 8839, 8849, 8861, 8863, 8867, 8887, 8893, 8923, 8929, 8933, 8941, 8951, 8963, 8969, 8971, 8999, 9001, 9007, 9011,
   9013, 9029, 9041, 9043, 9049, 9059, 9067, 9091, 9103, 9109, 9127, 9133, 9137, 9151, 9157, 9161, 9173, 9181, 9187, 9199,
   9203, 9209, 9221, 9227, 9239, 9241, 9257, 9277, 9281, 9283, 9293, 9311, 9319, 9323, 9337, 9341, 9343, 9349, 9371, 9377,
   9391, 9397, 9403, 9413, 9419, 9421, 9431, 9433, 9437, 9439, 9461, 9463, 9467, 9473, 9479, 9491, 9497, 9511, 9521, 9533,
   9539, 9547, 9551, 9587, 9601, 9613, 9619, 9623, 9629, 9631, 9643, 9649, 9661, 9677, 9679, 9689, 9697, 9719, 9721, 9733,
   9739, 9743, 9749, 9767, 9769, 9781, 9787, 9791, 9803, 9811, 9817, 9829, 9833, 9839, 9851, 9857, 9859, 9871, 9883, 9887,
   9901, 9907, 9923, 9929, 9931, 9941, 9949, 9967, 9973, 10007, 10009, 10037, 10039, 10061, 10067, 10069, 10079, 10091, 10093, 10099,
   10103, 10111, 10133, 10139, 10141, 10151, 10159, 10163, 10169, 10177, 10181, 10193, 10211, 10223, 10243, 10247, 10253, 10259, 10267, 10271,
   10273, 10289, 10301, 10303, 10313, 10321, 10331, 10333, 10337, 10343, 10357, 10369, 10391, 10399, 10427, 10429, 10433, 10453, 10457, 10459,
   10463, 10477, 10487, 10499, 10501, 10513, 10529, 10531, 10559, 10567, 10589, 10597, 10601, 10607, 10613, 10627, 10631, 10639, 10651, 10657,
   10663, 10667, 10687, 10691, 10709, 10711, 10723, 10729, 10733, 10739, 10753, 10771, 10781, 10789, 10799, 10831, 10837, 10847, 10853, 10859,
   10861, 10867, 10883, 10889, 10891, 10903, 10909, 10937, 10939, 10949, 10957, 10973, 10979, 10987, 10993, 11003, 11027, 11047, 11057, 11059,
   11069, 11071, 11083, 11087, 11093, 11113, 11117, 11119, 11131, 11149, 11159, 11161, 11171, 11173, 11177, 11197, 11213, 11239, 11243, 11251,
   11257, 11261, 11273, 11279, 11287, 11299, 11311, 11317, 11321, 11329, 11351, 11353, 11369, 11383, 11393, 11399, 11411, 11423, 11437, 11443,
   11447, 11467, 11471, 11483, 11489, 11491, 11497, 11503, 11519, 11527, 11549, 11551, 11579, 11587, 11593, 11597, 11617, 11621, 11633, 11657,
   11677, 11681, 11689, 11699, 11701, 11717, 11719, 11731, 11743, 11777, 11779, 11783, 11789, 11801, 11807, 11813, 11821, 11827, 11831, 11833,
   11839, 11863, 11867, 11887, 11897, 11903, 11909, 11923, 11927, 11933, 11939, 11941, 11953, 11959, 11969, 11971, 11981, 11987, 12007, 12011,
   12037, 12041, 12043, 12049, 12071, 12073, 12097, 12101, 12107, 12109, 12113, 12119, 12143, 12149, 12157, 12161, 12163, 12197, 12203, 12211,
   12227, 12239, 12241, 12251, 12253, 12263, 12269, 12277, 12281, 12289, 12301, 12323, 12329, 12343, 12347, 12373, 12377, 12379, 12391, 12401,
   12409, 12413, 12421, 12433, 12437, 12451, 12457, 12473, 12479, 12487, 12491, 12497, 12503, 12511, 12517, 12527, 12539, 12541, 12547, 12553,
   12569, 12577, 12583, 12589, 12601, 12611, 12613, 12619, 12637, 12641, 12647, 12653, 12659, 12671, 12689, 12697, 12703, 12713, 12721, 12739,
   12743, 12757, 12763, 12781, 12791, 12799, 12809, 12821, 12823, 12829, 12841, 12853, 12889, 12893, 12899, 12907, 12911, 12917, 12919, 12923,
   12941, 12953, 12959, 12967, 12973, 12979, 12983, 13001, 13003, 13007, 13009, 13033, 13037, 13043, 13049, 13063, 13093, 13099, 13103, 13109,
   13121, 13127, 13147, 13151, 13159, 13163, 13171, 13177, 13183, 13187, 13217, 13219, 13229, 13241, 13249, 13259, 13267, 13291, 13297, 13309,
   13313, 13327, 13331, 13337, 13339, 13367, 13381, 13397, 13399, 13411, 13417, 13421, 13441, 13451, 13457, 13463, 13469, 13477, 13487, 13499,
   13513, 13523, 13537, 13553, 13567, 13577, 13591, 13597, 13613, 13619, 13627, 13633, 13649, 13669, 13679, 13681, 13687, 13691, 13693, 13697,
   13709, 13711, 13721, 13723, 13729, 13751, 13757, 13759, 13763, 13781, 13789, 13799, 13807, 13829, 13831, 13841, 13859, 13873, 13877, 13879,
   13883, 13901, 13903, 13907, 13913, 13921, 13931, 13933, 13963, 13967, 13997, 13999, 14009, 14011, 14029, 14033, 14051, 14057, 14071, 14081,
   14083, 14087, 14107, 14143, 14149, 14153, 14159, 14173, 14177, 14197, 14207, 14221, 14243, 14249, 14251, 14281, 14293, 14303, 14321, 14323,
   14327, 14341, 14347, 14369, 14387, 14389, 14401, 14407, 14411, 14419, 14423, 14431, 14437, 14447, 14449, 14461, 14479, 14489, 14503, 14519,
   14533, 14537, 14543, 14549, 14551, 14557, 14561, 14563, 14591, 14593, 14621, 14627, 14629, 14633, 14639, 14653, 14657, 14669, 14683, 14699,
   14713, 14717, 14723, 14731, 14737, 14741, 14747, 14753, 14759, 14767, 14771, 14779, 14783, 14797, 14813, 14821, 14827, 14831, 14843, 14851,
   14867, 14869, 14879, 14887, 14891, 14897, 14923, 14929, 14939, 14947, 14951, 14957, 14969, 14983, 15013, 15017, 15031, 15053, 15061, 15073,
   15077, 15083, 15091, 15101, 15107, 15121, 15131, 15137, 15139, 15149, 15161, 15173, 15187, 15193, 15199, 15217, 15227, 15233, 15241, 15259,
   15263, 15269, 15271, 15277, 15287, 15289, 15299, 15307, 15313, 15319, 15329, 15331, 15349, 15359, 15361, 15373, 15377, 15383, 15391, 15401,
   15413, 15427, 15439, 15443, 15451, 15461, 15467, 15473, 15493, 15497, 15511, 15527, 15541, 15551, 15559, 15569, 15581, 15583, 15601, 15607,
   15619, 15629, 15641, 15643, 15647, 15649, 15661, 15667, 15671, 15679, 15683, 15727, 15731, 15733, 15737, 15739, 15749, 15761, 15767, 15773,
   15787, 15791, 15797, 15803, 15809, 15817, 15823, 15859, 15877, 15881, 15887, 15889, 15901, 15907, 15913, 15919, 15923, 15937, 15959, 15971,
   15973, 15991, 16001, 16007, 16033, 16057, 16061, 16063, 16067, 16069, 16073, 16087, 16091, 16097, 16103, 16111, 16127, 16139, 16141, 16183,
   16187, 16189, 16193, 16217, 16223, 16229, 16231, 16249, 16253, 16267, 16273, 16301, 16319, 16333, 16339, 16349, 16361, 16363, 16369, 16381,
   16411, 16417, 16421, 16427, 16433, 16447, 16451, 16453, 16477, 16481, 16487, 16493, 16519, 16529, 16547, 16553, 16561, 16567, 16573, 16603,
   16607, 16619, 16631, 16633, 16649, 16651, 16657, 16661, 16673, 16691, 16693, 16699, 16703, 16729, 16741, 16747, 16759, 16763, 16787, 16811,
   16823, 16829, 16831, 16843, 16871, 16879, 16883, 16889, 16901, 16903, 16921, 16927, 16931, 16937, 16943, 16963, 16979, 16981, 16987, 16993,
   17011, 17021, 17027, 17029, 17033, 17041, 17047, 17053, 17077, 17093, 17099, 17107, 17117, 17123, 17137, 17159, 17167, 17183, 17189, 17191,
   17203, 17207, 17209, 17231, 17239, 17257, 17291, 17293, 17299, 17317, 17321, 17327, 17333, 17341, 17351, 17359, 17377, 17383, 17387, 17389,
   17393, 17401, 17417, 17419, 17431, 17443, 17449, 17467, 17471, 17477, 17483, 17489, 17491, 17497, 17509, 17519, 17539, 17551, 17569, 17573,
   17579, 17581, 17597, 17599, 17609, 17623, 17627, 17657, 17659, 17669, 17681, 17683, 17707, 17713, 17729, 17737, 17747, 17749, 17761, 17783, 17789, 17791, 17807, 17827, 17837, 17839, 17851, 17863];
   
   #[test]
 fn check_sum(){
 let mut sum = 0u64;
 
   for i in PRIME_INV_64{
    sum = sum.wrapping_add(i)
   }
   assert_eq!(sum,14247534216423351310);
   sum = 0;
   
   for i in PRIME_INV_128{
     sum = sum.wrapping_add(i as u64).wrapping_add((i>>64) as u64)
   }
   assert_eq!(sum, 2833601251621510603);
   sum = 0;
   for i in PRIMELIST{
        sum = sum.wrapping_add(i as u64);
   } 
   assert_eq!(sum,17120309)
   
 }  
