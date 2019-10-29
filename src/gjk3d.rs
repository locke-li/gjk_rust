use crate::vector::Float3;
use crate::minkowski_sum_helper;

pub fn check(polya_:&[Float3], polyb_:&[Float3], ab_:Float3) -> (bool, Float3) {
    const max_iteration: i32 = 32;
    let mut iteration = 0;
    if ab_.is_zero() {
        ab_ = Float3 {x:1.0,y:1.0,z:1.0};
    }
    let s0 = minkowski_sum_helper::support(polya_, polyb_, &ab_);
    let s1 = minkowski_sum_helper::support(polya_, polyb_, -s0.v);

    (true, Float3 {x:0.0, y:0.0, z:0.0})
}