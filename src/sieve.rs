pub(crate) fn erasto_sieve(sup: usize) -> Vec<u32> {
    let ndxlmt = (sup - 3) / 2 + 1;
    let bfsz = ((sup - 3) / 2) / 32 + 1;

    let mut cmpsts = vec![0u32; bfsz];
    let sqrtndxlmt = ((sup as f64).sqrt() as usize - 3) / 2 + 1;

    for ndx in 0..sqrtndxlmt {
        if (cmpsts[ndx >> 5] & (1u32 << (ndx & 31))) == 0 {
            let p = ndx + ndx + 3;

            let mut cullpos = (p * p - 3) / 2;
            while cullpos < ndxlmt {
                unsafe {
                    let cptr = cmpsts.get_unchecked_mut(cullpos >> 5);
                    *cptr |= 1u32 << (cullpos & 31);
                }
                cullpos += p;
            }
        }
    }

    cmpsts
}

pub(crate) fn prime_count(sup: usize, data: &[u32]) -> u32 {
    if sup < 2 {
        return 0;
    }
    let mut corrector = 0u32;
    if sup & 1 == 1 {
        corrector = 1;
    }
    ((sup as u64 / 64u64) * 32u64 + (sup as u64 % 64) / 2
        - data.iter().map(|&x| x.count_ones() as u64).sum::<u64>()) as u32
        + corrector
}

pub(crate) fn prime_list_32(inf: usize, sup: usize, data: &[u32]) -> Vec<u32> {
    if sup < 2 {
        return vec![];
    }
    let mut primes = vec![2];
    let ndxlmt = (sup - 3) / 2 + 1;
    let lo: isize = (inf as isize - 3) / 2 + 1; // lo = -1;
    let temp = (lo..ndxlmt as isize)
        .into_iter()
        .filter_map(move |i| {
            if i < 0 {
                Some(2)
            } else if data[i as usize >> 5] & (1u32 << (i & 31)) == 0 {
                Some((i + i + 3) as u32)
            } else {
                None
            }
        })
        .collect::<Vec<u32>>();
    primes.extend_from_slice(&temp[..]);
    primes
}
