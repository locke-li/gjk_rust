use crate::vector::Float3;
use crate::minkowski_sum_helper;
type Point = minkowski_sum_helper::MinkowskiSumPoint<Float3>;

pub fn check(polya_:&[Float3], polyb_:&[Float3], ab_:&Float3) -> Option<(bool, Float3)> {
    const MAX_ITERATION: i32 = 32;
    let mut iteration = 0;
    let d = if ab_.is_zero() { &Float3 {x:1.0,y:1.0,z:1.0} } else { ab_ };
    let mut s0 = minkowski_sum_helper::support(polya_, polyb_, &d);
    let mut s1 = minkowski_sum_helper::support(polya_, polyb_, &-s0.v);
    let d = s1.v - s0.v;
    let d = Float3::triple_cross(&d, &-s0.v, &d);
    let mut s2 = minkowski_sum_helper::support(polya_, polyb_, &d);
    let mut d = Float3::triangle_normal(&s0.v, &s1.v, &s2.v);
    if d.dot(&s0.v) > 0.0 { let temp=s1; s1=s2; s2=temp; d=-d };
    let mut s3: Point;
    loop {
        s3 = minkowski_sum_helper::support(polya_, polyb_, &d);
        if s3.v.dot(&d) < 0.0 { break no_collision(polya_, polyb_, s0, s1, s2, s3, d) }
        let n0 = Float3::triangle_normal(&s0.v, &s1.v, &s3.v);
        let n1 = Float3::triangle_normal(&s2.v, &s0.v, &s3.v);
        let n2 = Float3::triangle_normal(&s1.v, &s2.v, &s3.v);
        let d0 = s3.v.dot(&n0);
        let d1 = s3.v.dot(&n1);
        let d2 = s3.v.dot(&n2);
        if d0 >= 0.0 && d1 >= 0.0 && d2 >= 0.0 {
            break epa(polya_, polyb_, s0, s1, s2, s3, d);
        }
        if d0 < 0.0 {
            s2 = s3;
            d = n0;
        }
        else if d1 < 0.0 {
            s1 = s3;
            d = n1;
        }
        else if d2 < 0.0 {
            s0 = s3;
            d = n2;
        }
        iteration += 1;
        if iteration >= MAX_ITERATION { break None }
    }
}

fn no_collision(
    polya_:&[Float3], polyb_:&[Float3], 
    mut s0:Point, mut s1:Point, mut s2:Point, mut s3:Point, mut d:Float3
    ) -> Option<(bool, Float3)> {

    Some((false, Float3::zero()))
}

fn epa(
    polya_:&[Float3], polyb_:&[Float3], 
    mut s0:Point, mut s1:Point, mut s2:Point, mut s3:Point, mut d:Float3
    ) -> Option<(bool, Float3)> {

    Some((true, Float3::zero()))
}