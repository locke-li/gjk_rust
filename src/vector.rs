#[derive(Copy, Clone)]
#[derive(PartialEq)]
#[derive(Debug)]
pub struct Float3 {
    pub x: f32,
    pub y: f32,
    pub z: f32,
}

impl Float3 {
    pub fn zero() -> Float3 {
        Float3 {x:0.0, y:0.0, z:0.0}
    }

    pub fn is_zero(&self) -> bool {
        self.x == 0.0 && self.y == 0.0
    }

    pub fn sqr_magnitude(&self) -> f32 {
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    pub fn magnitude(&self) -> f32 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    pub fn normalize(&mut self) {
        let m = self.magnitude();
        self.x /= m;
        self.y /= m;
        self.z /= m;
    }

    pub fn dot(&self, b:&Float3) -> f32 {
        self.x * b.x + self.y * b.y + self.z * b.z
    }

    pub fn cross(&self, b:&Float3) -> Float3 {
        Float3 {
            x: self.y * b.z - self.z * b.y,
            y: self.x * b.z - self.z * b.x,
            z: self.x * b.y - self.y * b.x,
        }
    }

    pub fn triangle_normal(&self, b:&Float3, c:&Float3) -> Float3 {
        Float3::cross(&(*b - *self), &(*c - *self))
    }

    pub fn triple_cross(&self, b:&Float3, c:&Float3) -> Float3 {
        Float3::cross(&Float3::cross(self, b), c)
    }

    pub fn scale(&mut self, p_:f32, b_:f32) {
        let f = p_ / b_;
        self.x *= f;
        self.y *= f;
        self.z *= f;
    }
}

use std::ops;

impl ops::Add for Float3 {
    type Output = Float3;

    fn add(self, rhs_: Float3) -> Float3 {
        Float3 {
            x: self.x + rhs_.x,
            y: self.y + rhs_.y,
            z: self.z + rhs_.z,
        }
    }
}

impl ops::Sub for Float3 {
    type Output = Float3;

    fn sub(self, rhs_: Float3) -> Float3 {
        Float3 {
            x: self.x - rhs_.x,
            y: self.y - rhs_.y,
            z: self.z - rhs_.z,
        }
    }
}

impl ops::Neg for Float3 {
    type Output = Float3;

    fn neg(self) -> Float3 {
        Float3 {x:-self.x, y:-self.y, z:-self.z}
    }
}

impl ops::Neg for &Float3 {
    type Output = Float3;

    fn neg(self) -> Float3 {
        Float3 {x:-self.x, y:-self.y, z:-self.z}
    }
}