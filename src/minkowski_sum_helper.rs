use crate::vector::Float3;

#[derive(Debug)]
pub struct MinkowskiSumPoint<T> {
    pub v: T,
    pub a: usize,
    pub b: usize,
}

impl<T> std::cmp::PartialEq for MinkowskiSumPoint<T> {
    fn eq(&self, rhs_:&Self) -> bool {
        self.a == rhs_.a && self.b == rhs_.b
    }
}

pub fn support(polya_:&[Float3], polyb_:&[Float3], d_:&Float3) -> MinkowskiSumPoint<Float3> {
    let a = support_in(polya_, d_);
    let b = support_in(polyb_, &-d_);
    println!("dir {} support {}", d_, polya_[a] - polyb_[b]);
    MinkowskiSumPoint {v: polya_[a] - polyb_[b], a:a, b:b}
}

fn support_in(poly_:&[Float3], d_:&Float3) -> usize {
    let mut ret = 0;
    let mut max: f32 = std::f32::MIN;
    for (i, v) in poly_.iter().enumerate() {
        let c = v.dot(d_);
        if c > max {
            max = c;
            ret = i;
        }
    }
    ret
}