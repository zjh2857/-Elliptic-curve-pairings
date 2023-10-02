use std::ops::{Add,Sub,Neg,Mul,Div};
use std::cmp::PartialEq;

const P : i128 = 60236498173;
const B : i128 = 7;
const BETA : i128 = 2;
const R : i128 = 60236253349;
const RINV : i128 = 44310963719;
const RBITS : i128 = 36;

#[derive(Copy, Clone,Debug)]
struct Fp<T> {
    x : T
}
impl Fp<i128> {
    // fn new(x : i128) -> Fp<i128> {
    //     Fp { x: x % P }
    // }
    fn powmod(&self, n : u64) -> Fp<i128> {
        let mut n = n;
        let mut x = self.x;
        let mut r = 1;
        while n > 0 {
            if n & 1 == 1 {
                r = (r * x) % P;
            }
            x = (x * x) % P;
            n >>= 1;
        }
        Fp { x: (r) }
    }
    fn powmodbits(&self, bits:Vec<u64>) -> Fp<i128> {
        let mut x = self.x;
        let mut r = 1;
        for i in 0..bits.len() {
            if bits[i] == 1 {
                r = (r * x) % R;
            }
            x = (x * x) % R;
        }
        Fp { x: (r) }
    }
}
impl Add for Fp<i128> {
    type Output = Fp<i128>;

    fn add(self, other: Fp<i128>) -> Fp<i128> {
        Fp { x: (self.x + other.x) % P }
    }
}

impl Sub for Fp<i128> {
    type Output = Fp<i128>;

    fn sub(self, other: Fp<i128>) -> Fp<i128> {
        Fp { x: (P + self.x - other.x) % P }
    }
}
    
impl Mul for Fp<i128> {
    type Output = Fp<i128>;

    fn mul(self, other: Fp<i128>) -> Fp<i128> {
        Fp { x: (self.x * other.x) % P }
    }
}


impl Neg for Fp<i128> {
    type Output = Fp<i128>;

    fn neg(self) -> Fp<i128> {
        Fp { x: P-self.x % P}
    }
}

impl PartialEq<Fp<i128>> for Fp<i128> {
    fn eq(&self, other: &Fp<i128>) -> bool {
        self.x == other.x
    }
}

impl Div for Fp<i128> {
    type Output = Fp<i128>;

    fn div(self, other: Fp<i128>) -> Fp<i128> {
        self * other.powmod(P as u64 - 2)
    }
}
    
#[derive(Copy, Clone,Debug)]
struct F12 {
    coff : [Fp<i128>; 12]
}


impl F12 {
    fn new() -> F12 {
        F12 { coff: [Fp { x: 0 }; 12] }
    }
    fn powmod(&self, n : u64) -> F12 {
        let mut n = n;
        let mut x = self.clone();
        let mut r = F12::from_i128(1);
        while n > 0 {
            if n & 1 == 1 {
                r = r * x;
            }
            x = x * x;
            n >>= 1;
        }
        r
    }
    fn powmodbits(&self, bits:Vec<u64>) -> F12 {
        let mut x = self.clone();
        let mut r = F12::from_i128(1);
        for i in 0..bits.len() {
            if bits[i] == 1 {
                r = r * x;
            }
            x = x * x;
        }
        r
    }
    fn from_fp(x : Fp<i128>) -> F12 {
        let mut r = F12::new();
        r.coff[0] = x;
        r
    }
    fn from_i128(x : i128) -> F12 {
        let mut r = F12::new();
        r.coff[0] = Fp { x: x };
        r
    }
    fn from_i128_vec(x : Vec<i128>) -> F12 {
        let mut r = F12::new();
        for i in 0..12 {
            r.coff[i] = Fp { x: x[i] };
        }
        r
    }
    fn gauss_elimination(&self, other : F12) -> F12{
        // Gauss
        let mut matrix:Vec<Vec<Fp<i128>>> = vec![vec![Fp{x:0};13];12];

        for i in 0..12 {
            for j in 0..12 {
                if i >= j {
                    matrix[i][j] = other.coff[i-j];
                } else {
                    matrix[i][j] = other.coff[i+12-j] * Fp { x: BETA };
                }
            }
        }
        for i in 0..12 {
            matrix[i][12] = self.coff[i];
        }

        let size = 12;
        for i in 0..size - 1 {
            if matrix[i][i].x == 0 {
                for j in i..size - 1 {
                    if matrix[j][i].x != 0 {
                        for k in 0..size + 1 {
                            let tmp = matrix[i][k];
                            matrix[i][k] = matrix[j][k];
                            matrix[j][k] = tmp;
                        }
                        break;
                    }
                }
            }
            for j in i..size - 1 {
                echelon(&mut matrix, i, j);
            }
        }

        for i in (1..size).rev() {
            eliminate(&mut matrix, i);
            // debug(&matrix);
        }
        let mut result: F12 = F12::new();
        for i in 0..size {
            result.coff[i] = matrix[i][size] / matrix[i][i];
        }
        result
    }
}

// fn debug(matrix : &Vec<Vec<Fp<i128>>>){
//     for i in 0..12 {
//         for j in 0..13 {
//             print!("{}", (matrix[i][j].x));
//         }
//         print!("\n");
//     }
//     print!("\n\n");
// }
fn echelon(matrix: &mut Vec<Vec<Fp<i128>>>, i: usize, j: usize) {
    let size = matrix.len();
    if matrix[i][i].x == 0 {
        
    } else {
        let factor = matrix[j + 1][i] / matrix[i][i];
        (i..size + 1).for_each(|k| {
            matrix[j + 1][k] = matrix[j + 1][k] - factor * matrix[i][k];
        });
    }
}

fn eliminate(matrix: &mut Vec<Vec<Fp<i128>>>, i: usize) {
    let size = matrix.len();
    if matrix[i][i].x == 0 {
    } else {
        for j in (1..i + 1).rev() {
            let factor = matrix[j - 1][i] / matrix[i][i];
            for k in (0..size + 1).rev() {
                matrix[j - 1][k] = matrix[j - 1][k] - factor * matrix[i][k];
            }
        }
    }
}
impl Add for F12 {
    type Output = F12;

    fn add(self, other: F12) -> F12 {
        let mut r = F12::new();
        for i in 0..12 {
            r.coff[i] = self.coff[i] + other.coff[i];
        }
        r
    }
}

impl Sub for F12 {
    type Output = F12;

    fn sub(self, other: F12) -> F12 {
        let mut r = F12::new();
        for i in 0..12 {
            r.coff[i] = self.coff[i] - other.coff[i];
        }
        r
    }
}

impl Mul for F12 {
    type Output = F12;

    fn mul(self, other: F12) -> F12 {
        let mut r = F12::new();
        for i in 0..12 {
            for j in 0..12 {
                if i+j < 12 {
                    r.coff[i+j] = r.coff[i+j] + self.coff[i] * other.coff[j];
                } else {
                    r.coff[i+j-12] = r.coff[i+j-12] + self.coff[i] * other.coff[j] * Fp { x: BETA };
                }
            }
        }
        r
    }
}

impl Neg for F12 {
    type Output = F12;

    fn neg(self) -> F12 {
        let mut r = F12::new();
        for i in 0..12 {
            r.coff[i] = -self.coff[i];
        }
        r
    }
}

impl Div for F12 {
    type Output = F12;

    fn div(self, other: F12) -> F12 {
        let r: F12 = self.gauss_elimination(other);
        r
    }
}

impl PartialEq<F12> for F12 {
    fn eq(&self, other: &F12) -> bool {
        for i in 0..12 {
            if self.coff[i] != other.coff[i] {
                return false;
            }
        }
        true
    }
}

    
#[allow(dead_code)]
#[derive(Copy, Clone,Debug)]
struct Point {
    x : Fp<i128>,
    y : Fp<i128>,
    z : Fp<i128>   
}

impl Point {
    fn scale(self, other: i128) -> Point {
        let mut r = Point { x: Fp { x: 0 }, y: Fp { x: 1 }, z: Fp { x: 0 } };
        let mut p = self;
        let mut q: i128 = other;
        while q > 0 {
            if q & 1 == 1 {
                r = r + p;
            }
            p = p.clone() + p.clone();
            q >>= 1;
        }
        r
    }    
}

impl Add for Point {
    type Output = Point;

    fn add(self, other: Point) -> Point {
        if self.z == (Fp { x: 0 }) {
            return other;
        }
        else if other.z == (Fp { x: 0 }) {
            return self;
        }
        else if self.x == other.x && self.y == -other.y {
            Point { x: Fp { x: 0 }, y: Fp { x: 1 }, z: Fp { x: 0 } }
        } 
        else if self.x == other.x && self.y == other.y {
            let s = (self.x * self.x * Fp { x: 3 }) / (self.y * Fp { x: 2 });
            let x = s * s - self.x - other.x;
            let y = s * (self.x - x) - self.y;

            Point { x: x, y: y, z: Fp { x: 1 } }
        } else {
            let s = (self.y - other.y) / (self.x - other.x);
            let x = s * s - self.x - other.x;
            let y = s * (self.x - x) - self.y;
            Point { x: x, y: y, z: Fp { x: 1 } }
        }
    }
}

impl Sub for Point {
    type Output = Point;

    fn sub(self, other: Point) -> Point {
        self + (-other)
    }
}

impl Neg for Point {
    type Output = Point;

    fn neg(self) -> Point {
        Point { x: self.x, y: -self.y, z: self.z }
    }
}

#[allow(dead_code)]
#[derive(Copy, Clone,Debug)]
struct PointExt {
    x : F12,
    y : F12,
    z : F12   
}

impl PointExt {
    // fn new(x : Fp<i128>, y : Fp<i128>) -> PointExt {
    //     PointExt { x: x, y: y, z: Fp { x: 1 } }
    // }
    // fn is_on_curve(&self) -> bool {
    //     self.y * self.y == self.x * self.x * self.x + Fp { x: B }
    // }
    fn scale(self, other: i128) -> PointExt {
        let mut r = PointExt { x: F12::new(), y: F12::new(), z: F12::new() };
        let mut p = self;
        let mut q: i128 = other;
        while q > 0 {
            if q & 1 == 1 {
                r = r + p;
            }
            p = p.clone() + p.clone();
            q >>= 1;
        }
        r
    }

}

impl Add for PointExt {
    type Output = PointExt;

    fn add(self, other: PointExt) -> PointExt {
        if self.z == F12::new() {
            return other;
        }
        else if other.z == F12::new() {
            return self;
        }
        else if self.x == other.x && self.y == -other.y {
            PointExt { x: F12::new(), y: F12::from_i128(1), z: F12::new() }
        } 
        else if self.x == other.x && self.y == other.y {
            let s = (self.x * self.x * F12::from_i128(3) ) / (self.y * F12::from_i128(2));
            let x = s * s - self.x - other.x;
            let y = s * (self.x - x) - self.y;
            PointExt { x: x, y: y, z: F12::from_i128(1) }
        } else {
            let s = (self.y - other.y) / (self.x - other.x);
            let x = s * s - self.x - other.x;
            let y = s * (self.x - x) - self.y;
            PointExt { x: x, y: y, z: F12::from_i128(1) }
        }
    }    
}

impl Neg for PointExt {
    type Output = PointExt;

    fn neg(self) -> PointExt {
        PointExt { x: self.x, y: -self.y, z: self.z }
    }
}

impl Sub for PointExt {
    type Output = PointExt;

    fn sub(self, other: PointExt) -> PointExt {
        self + (-other)
    }
}
    
fn miller_loop(p: PointExt, q: PointExt) -> F12 {
    let mut f: F12 = F12::from_i128(1);
    let mut mp = p;
    // let mut g: F12 = F12::from_i128(1);
    // let mut mq = q;


    let bits = vec![1,0,1,0,0,1,0,1,0,0,0,1,0,0,1,0,0,0,1,1,1,0,1,0,0,1,1,0,0,0,0,0,0,1,1,1];

    for i in 1..RBITS {
        let lmp2 = mp.x * mp.x * F12::from_i128(3) / (mp.y * F12::from_i128(2)) * (q.x - mp.x) + mp.y - q.y;
        let l2mp =  q.x - mp.x;
        f = f * f * lmp2 / l2mp;
        mp = mp + mp;
        
        if bits[i as usize] == 1 {
            let s = (mp.y - p.y) / (mp.x - p.x);
            let lmpp = s * (q.x - p.x) + p.y - q.y;
            let lmp = q.x - mp.x;
            mp = p + mp;
            f = f * lmpp / lmp;
        }
    }

    // for i in 1..RBITS {
    //     let lmq2 = mq.x * mq.x * F12::from_i128(3) / (mq.y * F12::from_i128(2)) * (p.x - mq.x) + mq.y;
    //     let l2mq =  p.x - mq.x;
    //     g = g * g * lmq2 / l2mq;
    //     mq = mq + mq;
        
    //     if bits[i as usize] == 1 {
    //         let s = (mq.y - q.y) / (mq.x - q.x);
    //         let lmqq = s * (p.x - q.x) + q.y ;
    //         let lmq = p.x - mq.x;
    //         mq = q + mq;
    //         g = g * lmqq / lmq;
    //     }
    // }
    let k12dr : Vec<u64> = [0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1].to_vec();

    f.powmodbits(k12dr)
}
#[cfg(test)]
mod tests {
    use ndarray_linalg::assert;

    use super::*;

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
        print!("{}\n", c.coff[0].x);
        let d: F12 = c / b;
        print!("{}\n", d.coff[0].x);
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
        let xa = F12::from_i128_vec(vec![59649468771, 27477382595, 54970258839, 24464396194, 50445504362, 6691406414, 15613152500, 57679389230, 36873964818, 49047050504, 37627375104,25144576512]);
        let ya = F12::from_i128_vec(vec![21585877989, 22868212910, 15172448012, 14052413555, 45604206563, 44997454193, 9764794103, 21639994320, 24355274013, 55237684856, 21354386475,56640594939]);
        let xa = F12::from_i128(49015138457);
        let ya = F12::from_i128(36336695411);
        let a = PointExt{x : xa, y : ya, z : F12::from_i128(1)};
        let xb = F12::from_i128_vec(vec![52649681422, 58330021724, 18361934063, 37544484019, 48378959847, 34344081429, 20873394922, 26784166228, 25856414820, 31733507710, 1579608448,25687127276]);
        let yb = F12::from_i128_vec(vec![36962872648, 28863545577, 49052618431, 14807993877, 27927561526, 25370717733, 16963613074, 45537344140, 11799967561, 30748944721, 10754770505,42033372585]);
        let b = PointExt{x : xb, y : yb, z : F12::from_i128(1)};
        let e = b + b;
        let c = miller_loop(a,b);
        let d = miller_loop(a,e);
        print!("{}\n", c.powmod(R as u64).coff[0].x);
        // print!("{}\n", (c*c).coff[0].x);
        print!("{}\n",(c * c / d).coff[0].x );
        assert!((c * c / d).coff[0].x == 1);
    }
}
