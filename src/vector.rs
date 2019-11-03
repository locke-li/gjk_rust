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
            y: self.z * b.x - self.x * b.z,
            z: self.x * b.y - self.y * b.x,
        }
    }

    pub fn triangle_normal(&self, b:&Float3, c:&Float3) -> Float3 {
        let l0 = *b - *self;
        let l1 = *c - *self;
        println!("{} {} {} => {} {}", self, b, c, l0, l1);
        Float3::cross(&l0, &l1)
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

impl std::fmt::Display for Float3 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "({},{},{})", self.x, self.y, self.z)
    }
}

impl std::ops::Add for Float3 {
    type Output = Float3;

    fn add(self, rhs_: Float3) -> Float3 {
        Float3 {
            x: self.x + rhs_.x,
            y: self.y + rhs_.y,
            z: self.z + rhs_.z,
        }
    }
}

impl std::ops::Sub for Float3 {
    type Output = Float3;

    fn sub(self, rhs_: Float3) -> Float3 {
        Float3 {
            x: self.x - rhs_.x,
            y: self.y - rhs_.y,
            z: self.z - rhs_.z,
        }
    }
}

impl std::ops::Neg for Float3 {
    type Output = Float3;

    fn neg(self) -> Float3 {
        Float3 {x:-self.x, y:-self.y, z:-self.z}
    }
}

impl std::ops::Neg for &Float3 {
    type Output = Float3;

    fn neg(self) -> Float3 {
        Float3 {x:-self.x, y:-self.y, z:-self.z}
    }
}