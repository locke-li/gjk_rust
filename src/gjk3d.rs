use crate::vector::Float3;
use crate::simplex_based_cd_helper::*;
type Point = MinkowskiSumPoint<Float3>;
type Frame = Frame3Simplex<Float3>;

const MAX_ITERATION: i32 = 32;

pub fn check(polya_:&[Float3], polyb_:&[Float3], ab_:&Float3, f_:&mut Frame) -> Result<bool, Error> {
    if polya_.len() < 4 || polyb_.len() < 4 { return Err(Error::InvalidInput) }
    let mut iteration = 0;
    let d = if ab_.is_zero() { &Float3 {x:1.0,y:1.0,z:1.0} } else { ab_ };
    let mut s0 = support(polya_, polyb_, &d);
    let mut s1 = support(polya_, polyb_, &-s0.v);
    let d = s1.v - s0.v;
    let d = plane_normal(&d, &-s0.v);
    let mut s2 = support(polya_, polyb_, &d);
    let mut d = Float3::triangle_normal(&s0.v, &s1.v, &s2.v);
    if d.dot(&s0.v) > 0.0 {
        std::mem::swap(&mut s1, &mut s2);
        d=-d;
    }
    let mut s3: Point;
    loop {
        s3 = support(polya_, polyb_, &d);
        println!("{}\n{} {} {} {}\n{}", d, s0, s1, s2, s3, s3.v.dot(&d));
        if s3.v.dot(&d) < 0.0 { break no_collision(polya_, polyb_, f_, s0, s1, s2, s3, d) }
        let n0 = Float3::triangle_normal(&s0.v, &s1.v, &s3.v);
        let n1 = Float3::triangle_normal(&s2.v, &s0.v, &s3.v);
        let n2 = Float3::triangle_normal(&s1.v, &s2.v, &s3.v);
        let d0 = s3.v.dot(&n0);
        let d1 = s3.v.dot(&n1);
        let d2 = s3.v.dot(&n2);
        println!("{} {} {} {} {} {}", n0, n1, n2, d0, d1, d2);
        if d0 >= 0.0 && d1 >= 0.0 && d2 >= 0.0 {
            break epa(polya_, polyb_, f_, s0, s1, s2, s3, d);
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
        if iteration >= MAX_ITERATION { break Err(Error::SimplexSearch) }
    }
}

fn no_collision(
    polya_:&[Float3], polyb_:&[Float3], f_:&mut Frame,
    mut s0:Point, mut s1:Point, mut s2:Point, mut s3:Point, mut d:Float3
) -> Result<bool, Error> {
    let mut iteration = 0;
    while d.dot(&(s3.v - s0.v)) > 0.0 {
        let n0 = Float3::triangle_normal(&s0.v, &s1.v, &s3.v);
        let n1 = Float3::triangle_normal(&s3.v, &s1.v, &s2.v);
        let n2 = Float3::triangle_normal(&s0.v, &s3.v, &s2.v);
        let d0 = n0.dot(&s3.v);
        let d0 = if d0 < 0.0 { 0.0 } else { d0 };
        let d1 = n1.dot(&s3.v);
        let d1 = if d1 < 0.0 { 0.0 } else { d1 };
        let d2 = n2.dot(&s3.v);
        let d2 = if d2 < 0.0 { 0.0 } else { d2 };
        let d0 = -(d0 * d0) / n0.sqr_magnitude();
        let d1 = -(d1 * d1) / n1.sqr_magnitude();
        let d2 = -(d2 * d2) / n2.sqr_magnitude();
        if d0 <= d1 && d0 <= d2 {
            s2 = s3;
        }
        else if d1 <= d0 && d1 <= d2 {
            s0 = s3;
        }
        else if d2 <= d0 && d2 <= d1 {
            s1 = s3;
        }
        else {
            println!("unexpected {} {} {}", d0, d1, d2);
        }
        d = Float3::triangle_normal(&s0.v, &s1.v, &s2.v);
        s3 = support(polya_, polyb_, &d);
        iteration += 1;
        if iteration >= MAX_ITERATION { return Err(Error::NearestSimplexSearch) }
    }
    calculate_mtv_from_nearest_feature(polya_, polyb_, f_, s0, s1, s2, &d, false)
}

fn epa(
    polya_:&[Float3], polyb_:&[Float3], f_:&mut Frame,
    mut s0:Point, mut s1:Point, mut s2:Point, mut s3:Point, mut d:Float3
) -> Result<bool, Error> {
    f_.cache_epa.clear();
    f_.cache_epa.push(EPA2Simplex::new(s0, s2, s1));
    Ok(true)
}