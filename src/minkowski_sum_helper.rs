use crate::vector::Float3;

pub struct MinkowskiSumPoint<T> {
    pub v: T,
    pub a: i32,
    pub b: i32,
}

impl<T> std::cmp::PartialEq for MinkowskiSumPoint<T> {
    fn eq(&self, rhs_:&Self) -> bool {
        self.a == rhs_.a && self.b == rhs_.b
    }
}

pub fn support(polya_:&[Float3], polyb_:&[Float3], d_:&Float3) -> MinkowskiSumPoint<Float3> {
    

    MinkowskiSumPoint {v: Float3{x:0.0, y:0.0, z:0.0}, a:0, b:0}
}