pub mod field;
mod parma;
pub mod point;


pub mod pairing;

#[cfg(test)]
mod tests {
    use crate::field::{F12,Fp};
    use crate::point::{Point,PointExt};
    use crate::pairing::miller_loop_t;
    #[test]
    fn it_works() {
        let a:Fp<i128> = Fp { x: (123) };
        let b:Fp<i128> = Fp { x: (456) };
        let c:Fp<i128> = a + b;
        assert!(c.x == 123 + 456)
    }
    #[test]
    fn test_point() {
        let a = Point{x : Fp { x: 50884343186 }, y : Fp { x: 7946871878 }, z : Fp { x: 1 }};
        let b = Point{x : Fp { x: 6997785916 }, y : Fp { x: 3993585893 }, z : Fp { x: 1 }};
        let c = a + b;
        assert!(c.x.x == 15752766690);
        let a = Point{x : Fp { x: 37548606025 }, y : Fp { x: 28221163070 }, z : Fp { x: 1 }};
        let b = Point{x : Fp { x: 18771124967 }, y : Fp { x: 46507760661 }, z : Fp { x: 1 }};
        let c = a.scale(2);
        assert!(c.x.x == b.x.x);
    }
    #[test]
    fn test_pointext() {
        let xa = F12::from_i128_vec(vec![59649468771, 27477382595, 54970258839, 24464396194, 50445504362, 6691406414, 15613152500, 57679389230, 36873964818, 49047050504, 37627375104,25144576512]);
        let ya = F12::from_i128_vec(vec![21585877989, 22868212910, 15172448012, 14052413555, 45604206563, 44997454193, 9764794103, 21639994320, 24355274013, 55237684856, 21354386475,56640594939]);
        let a = PointExt{x : xa, y : ya, z : F12::from_i128(1)};
        let xb = F12::from_i128_vec(vec![1149259122, 43083045335, 36827928929, 40075706874, 4017825568, 20662801182, 6573426867, 58613218293, 25530413846, 21940044173, 44570367065,44874346429]);
        let yb = F12::from_i128_vec(vec![10754698008, 41486990217, 6932221709, 10941858426, 43830981973, 50886118260, 56425172390, 39396282642, 19750227184, 45902016622, 55187846690,25299303283]);
        // let scale = 3;
        let b = PointExt{x : xb, y : yb, z : F12::from_i128(1)};
        let c = a.scale(3);

        assert!(c.x.coff[0].x == b.x.coff[0].x);
        let d = a.scale(4) - a;
        assert!(d.x.coff[0].x == c.x.coff[0].x);
    }   
    #[test]
    fn test_f12() {
        let a: F12 = F12 { coff: [Fp { x: 1 }, Fp { x: 2 }, Fp { x: 3 }, Fp { x: 4 }, Fp { x: 5 }, Fp { x: 6 }, Fp { x: 7 }, Fp { x: 8 }, Fp { x: 9 }, Fp { x: 10 }, Fp { x: 11 }, Fp { x: 12 }] };
        let b: F12 = F12 { coff: [Fp { x: 1919 }, Fp { x: 2 }, Fp { x: 3 }, Fp { x: 810 }, Fp { x: 5 }, Fp { x: 6 }, Fp { x: 7 }, Fp { x: 8 }, Fp { x: 9 }, Fp { x: 10 }, Fp { x: 11 }, Fp { x: 12 }] };
        let c: F12 = a * b;
        let d: F12 = c / b;
        assert!(d.coff[0].x == 1);
        assert!(d.coff[1].x == 2);
        assert!(d.coff[2].x == 3);
        assert!(d.coff[3].x == 4);
        assert!(d.coff[4].x == 5);
        assert!(d.coff[5].x == 6);
        assert!(d.coff[6].x == 7);
        assert!(d.coff[7].x == 8);
        assert!(d.coff[8].x == 9);
    } 
    #[test]
    fn test_miller_loop(){
        let xa = F12::from_i128(49015138457);
        let ya = F12::from_i128(36336695411);
        let a = PointExt{x : xa, y : ya, z : F12::from_i128(1)};
        let xb = F12::from_i128_vec(vec![52649681422, 58330021724, 18361934063, 37544484019, 48378959847, 34344081429, 20873394922, 26784166228, 25856414820, 31733507710, 1579608448,25687127276]);
        let yb = F12::from_i128_vec(vec![36962872648, 28863545577, 49052618431, 14807993877, 27927561526, 25370717733, 16963613074, 45537344140, 11799967561, 30748944721, 10754770505,42033372585]);
        let b = PointExt{x : xb, y : yb, z : F12::from_i128(1)};
        let e = b.scale(3);
        let c = miller_loop_t(a,b);
        let d = miller_loop_t(a,e);
        assert!((c.powmod(3) / d).coff[0].x == 1);
    }
}
