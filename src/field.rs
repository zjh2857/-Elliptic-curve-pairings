use std::ops::{Add,Sub,Neg,Mul,Div};

use crate::parma::{P,BETA};


#[derive(Copy, Clone,Debug)]
pub struct Fp<T> {
    pub x : T
}

impl Fp<i128> {
    // fn new(x : i128) -> Fp<i128> {
    //     Fp { x: x % P }
    // }
    pub fn powmod(&self, n : u64) -> Fp<i128> {
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
    // pub fn powmodbits(&self, bits:Vec<u64>) -> Fp<i128> {
    //     let mut x = self.x;
    //     let mut r = 1;
    //     for i in 0..bits.len() {
    //         if bits[i] == 1 {
    //             r = (r * x) % R;
    //         }
    //         x = (x * x) % R;
    //     }
    //     Fp { x: (r) }
    // }
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
pub struct F12 {
    pub coff : [Fp<i128>; 12]
}

impl F12 {
    pub fn new() -> F12 {
        F12 { coff: [Fp { x: 0 }; 12] }
    }
    pub fn powmod(&self, n : u64) -> F12 {
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
    pub fn powmodbits(&self, bits:Vec<u64>) -> F12 {
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
    pub fn from_fp(x : Fp<i128>) -> F12 {
        let mut r = F12::new();
        r.coff[0] = x;
        r
    }
    pub fn from_i128(x : i128) -> F12 {
        let mut r = F12::new();
        r.coff[0] = Fp { x: x % P};
        r
    }
    pub fn from_i128_vec(x : Vec<i128>) -> F12 {
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
        if other == F12::new() {
            panic!("Divide by zero");
        }
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


