use number_theory::NumberTheory;

#[test]
fn sq_rt() {
    assert_eq!(0u8.sqrt().0, 0); // All n < 256 truncate to u8
    assert_eq!(1u8.sqrt().0, 1); //

    for i in 0..16u8 {
        let sqr = i * i;
        let rt = sqr.sqrt().0;
        assert_eq!(rt, i);
        assert_eq!(rt * rt, sqr);
    }

    for i in 16u16..256u16 {
        let sqr = i * i;
        let rt = sqr.sqrt().0;
        assert_eq!(rt, i);
        assert_eq!(rt * rt, sqr);
    }

    for i in 256u32..65536u32 {
        let sqr = i * i;
        let rt = sqr.sqrt().0;
        assert_eq!(rt, i);
        assert_eq!(rt * rt, sqr);
    }

    //assert_eq!
}
