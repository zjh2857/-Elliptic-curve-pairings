use crate::field::{F12,Fp};
use std::ops::{Add,Sub,Neg};
#[allow(dead_code)]
#[derive(Copy, Clone,Debug)]
pub struct Point {
    pub x : Fp<i128>,
    pub y : Fp<i128>,
    pub z : Fp<i128>   
}

impl Point {
    pub fn scale(self, other: i128) -> Point {
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
pub struct PointExt {
    pub x : F12,
    pub y : F12,
    pub z : F12   
}

impl PointExt {
    // fn new(x : Fp<i128>, y : Fp<i128>) -> PointExt {
    //     PointExt { x: x, y: y, z: Fp { x: 1 } }
    // }
    // fn is_on_curve(&self) -> bool {
    //     self.y * self.y == self.x * self.x * self.x + Fp { x: B }
    // }
    pub fn scale(self, other: i128) -> PointExt {
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

