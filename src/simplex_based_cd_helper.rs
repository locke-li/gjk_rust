use crate::vector::Float3;
type Point<T> = MinkowskiSumPoint<T>;

#[derive(Debug)]
pub enum Error {
    InvalidInput,
    SimplexSearch,
    NearestSimplexSearch,
    NearestFeatureSearch,
    EPA,
}

// region MinkowskiSumPoint

#[derive(Debug, Clone)]
pub struct MinkowskiSumPoint<T> {
    pub v: T,
    pub a: usize,
    pub b: usize,
}

impl<T> std::fmt::Display for MinkowskiSumPoint<T> where T: std::fmt::Display {
    fn fmt (&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} ({}-{})", self.v, self.a, self.b)
    }
}

impl<T> std::cmp::PartialEq for MinkowskiSumPoint<T> {
    fn eq(&self, rhs_:&Self) -> bool {
        self.a == rhs_.a && self.b == rhs_.b
    }
}

impl<T> Default for MinkowskiSumPoint<T> where T: Default {
    fn default() -> Self { MinkowskiSumPoint {v:Default::default(), a:0, b:0}}
}

// endregion

// region Frame3Simplex

#[derive(Debug)]
pub struct Frame3Simplex<T> {
    pub mtv: T,
    pub min_dist: f32,
    pub s0: Point<T>,
    pub s1: Point<T>,
    pub s2: Point<T>,
    pub closest_a: T,
    pub closest_b: T,
    pub cache_epa: Vec<EPA2Simplex<T>>,
    pub candidate_a: Vec<usize>,
    pub candidate_b: Vec<usize>,
}

impl<T> Frame3Simplex<T> where T : Default {
    pub fn new() -> Frame3Simplex<T> {
        Frame3Simplex {
            mtv: Default::default(),
            min_dist: 0.0,
            s0: Default::default(),
            s1: Default::default(),
            s2: Default::default(),
            closest_a: Default::default(),
            closest_b: Default::default(),
            cache_epa: Vec::new(),
            candidate_a: Vec::new(),
            candidate_b: Vec::new(),
        }
    }
}

impl<T> Frame3Simplex<T> {
    pub fn cache(&mut self, v0:Point<T>, v1:Point<T>, v2:Point<T>) {
        self.s0 = v0;
        self.s1 = v1;
        self.s2 = v2;
    }
}

// endregion

// region EPA2Simplex

#[derive(Debug, Clone)]
pub struct EPA2Simplex<T> {
    pub v0: Point<T>,
    pub v1: Point<T>,
    pub v2: Point<T>,
    pub n: T,
    pub p: T,
    pub e: T,
    pub d0: f32,
    pub d1: f32,
    pub d2: f32,
    pub d_sqr: f32,
}

impl EPA2Simplex<Float3> {
    pub fn new(v0_:&Point<Float3>, v1_:&Point<Float3>, v2_:&Point<Float3>) -> EPA2Simplex<Float3> {
        let n = Float3::triangle_normal(&v0_.v, &v1_.v, &v2_.v);
        let e0 = v1_.v - v0_.v;
        let n0 = Float3::cross(&e0, &n);
        let e1 = v2_.v - v1_.v;
        let n1 = Float3::cross(&e1, &n);
        let e2 = v0_.v - v2_.v;
        let n2 = Float3::cross(&e2, &n);
        let mut d0 = -n0.dot(&v0_.v);
        let mut d1 = -n1.dot(&v1_.v);
        let mut d2 = -n2.dot(&v2_.v);
        let e: Float3;
        let mut p: Float3;
        if d0 >= 0.0 && !e0.is_zero() {
            d2 = 0.0;
            e = e0;
            d0 = e.dot(&v0_.v);
            d1 = e0.sqr_magnitude();
            p = v0_.v + e0.scaled(d0, d1);
        }
        else if d1 >= 0.0 && !e1.is_zero() {
            d2 = 1.0;
            e = e1;
            d0 = e.dot(&v1_.v);
            d1 = e1.sqr_magnitude();
            p = v1_.v + e1.scaled(d0, d1);
        }
        else if d2 >= 0.0 && !e2.is_zero() {
            d2 = 2.0;
            e = e2;
            d0 = e.dot(&v2_.v);
            d1 = e2.sqr_magnitude();
            p = v2_.v + e2.scaled(d0, d1);
        }
        else {
            d2 = -1.0;
            e = e0;
            p = n;
            p.scale(n.dot(&v0_.v), n.sqr_magnitude());
        }
        EPA2Simplex {
            v0: v0_.clone(),
            v1: v1_.clone(),
            v2: v2_.clone(),
            n: n,
            p: p,
            e: e,
            d0: d0,
            d1: d1,
            d2: d2,
            d_sqr: p.sqr_magnitude(),
        }
    }
}

// endregion

pub fn plane_normal(v0_:&Float3, v1_:&Float3) -> Float3 {
    let cross01 = Float3::cross(v0_, v1_);
    if cross01.is_zero() {
        let x = v0_.x.abs();
        let y = v0_.y.abs();
        let z = v0_.z.abs();
        Float3::triple_cross(
            v0_,
            &if x < y && x < z { Float3 {x:1.0,y:0.0,z:0.0} } else
            if y < x && y < z { Float3 {x:0.0,y:1.0,z:0.0} } else
            { Float3 {x:0.0,y:0.0,z:1.0} },
            v0_,
        )
    }
    else {
        Float3::cross(&cross01, v0_)
    }
}

pub fn support(polya_:&[Float3], polyb_:&[Float3], d_:&Float3, f_:&mut Frame3Simplex<Float3>) -> MinkowskiSumPoint<Float3> {
    support_in(polya_, d_, &mut f_.candidate_a);
    support_in(polyb_, &-d_, &mut f_.candidate_b);
    let mut min_a = 0;
    let mut min_b = 0;
    let mut v = Default::default();
    let mut min = std::f32::MAX;
    for a in &f_.candidate_a {
        for b in &f_.candidate_b {
            let c = polya_[*a] - polyb_[*b];
            let n = c.sqr_magnitude();
            if n < min {
                min = n;
                min_a = *a;
                min_b = *b;
                v = c;
            }
        }
    }
    println!("dir {} support {}", d_, v);
    MinkowskiSumPoint {v:v, a:min_a, b:min_b}
}

fn support_in(poly_:&[Float3], d_:&Float3, candidate_:&mut Vec<usize>) {
    candidate_.clear();
    let mut max: f32 = std::f32::MIN;
    for (i, v) in poly_.iter().enumerate() {
        let c = v.dot(d_);
        if c > max {
            max = c;
            candidate_.clear();
        }
        if c >= max {
            candidate_.push(i);
        }
    }
}

fn better_support_in_candidate(polya_:&[Float3], polyb_:&[Float3], n_:&Float3,
    s1_:&mut Point<Float3>, s2_:&mut Point<Float3>, f_:&Frame3Simplex<Float3>
) -> bool {
    if f_.candidate_a.len() * f_.candidate_b.len() < 3 { return false; }
    let mut s3: Point<Float3> = Default::default();
    let mut max = 0.0;
    for a in &f_.candidate_a {
        for b in &f_.candidate_b {
            let v = polya_[*a] - polyb_[*b];
            let n = n_.dot(&(v - s1_.v));
            if n > max {
                max = n;
                s3.a = *a;
                s3.b = *b;
                s3.v = v;
            }
        }
    }
    if max > 0.0 {
        std::mem::swap(s1_, s2_);
        *s1_ = s3;
    }
    max > 0.0
}

pub fn calculate_mtv_from_nearest_feature(polya_:&[Float3], polyb_:&[Float3], f_:&mut Frame3Simplex<Float3>, 
    mut s0:Point<Float3>, mut s1:Point<Float3>, mut s2:Point<Float3>, d:&Float3
) -> Result<bool, Error> {
    let mut iteration = 0;
    const MAX_ITERATION: i32 = 16;
    loop {
        let n0 = d.cross(&(s1.v-s0.v));
        let n1 = d.cross(&(s2.v-s1.v));
        let n2 = d.cross(&(s0.v-s2.v));
        let d0 = -n0.dot(&s0.v);
        let d1 = -n1.dot(&s1.v);
        let d2 = -n2.dot(&s2.v);
        if d0 < 0.0 && d1 < 0.0 && d2 < 0.0 {
            f_.cache(s0, s1, s2);
            f_.mtv_from_face_case(&d, polya_, polyb_);
            break Ok(false)
        }
        else if d0 >= 0.0 && !n0.is_zero() && (d0 == 0.0 || !better_support_in_candidate(polya_, polyb_, &n0, &mut s1, &mut s2, f_)) {
            f_.mtv_from_edge_case(&s0, &s1, polya_, polyb_);
            f_.cache(s0, s1, s2);
            break Ok(false)
        }
        else if d1 >= 0.0 && !n1.is_zero() && (d1 == 0.0 || !better_support_in_candidate(polya_, polyb_, &n1, &mut s0, &mut s2, f_)) {
            f_.mtv_from_edge_case(&s1, &s2, polya_, polyb_);
            f_.cache(s0, s1, s2);
            break Ok(false)
        }
        else if d2 >= 0.0 && !n2.is_zero() && (d2 == 0.0 || !better_support_in_candidate(polya_, polyb_, &n0, &mut s0, &mut s1, f_)) {
            f_.mtv_from_edge_case(&s2, &s0, polya_, polyb_);
            f_.cache(s0, s1, s2);
            break Ok(false)
        }
        iteration += 1;
        if iteration >= MAX_ITERATION { break Err(Error::NearestFeatureSearch) }
    }
}

impl Frame3Simplex<Float3> {
    fn normalize(&mut self) {
        self.min_dist = self.mtv.magnitude();
        if self.min_dist > 0.0 { self.mtv.scale(1.0, self.min_dist); }
    }

    pub fn mtv_from_epa(&mut self, e_:&EPA2Simplex<Float3>, polya_:&[Float3], polyb_:&[Float3]) {
        self.s0 = e_.v0.clone();
        self.s1 = e_.v1.clone();
        self.s2 = e_.v2.clone();
        match e_.d2 as i32 {
            0 => self.mtv_from_edge_case_precomputed(polya_, polyb_, &e_.v0, &e_.v1, e_.d0, e_.d1, &e_.n),
            1 => self.mtv_from_edge_case_precomputed(polya_, polyb_, &e_.v1, &e_.v2, e_.d0, e_.d1, &e_.n),
            2 => self.mtv_from_edge_case_precomputed(polya_, polyb_, &e_.v2, &e_.v0, e_.d0, e_.d1, &e_.n),
            _ => self.mtv_from_face_case(&e_.p, polya_, polyb_),
        }
    }

    fn mtv_from_edge_case_precomputed(&mut self, polya_:&[Float3], polyb_:&[Float3],
        s0_:&Point<Float3>, s1_:&Point<Float3>, d0_:f32, d1_:f32, mtv_:&Float3
    ) {
        self.closest_a = Float3::lerp_clamp(&polya_[s0_.a], &polya_[s1_.a], d0_, d1_);
        self.closest_b = Float3::lerp_clamp(&polyb_[s0_.b], &polyb_[s1_.b], d0_, d1_);
        self.mtv = *mtv_;
        self.normalize();
    }

    fn mtv_from_edge_case(&mut self, s0_:&Point<Float3>, s1_:&Point<Float3>, polya_:&[Float3], polyb_:&[Float3]) {
        let e = s1_.v - s0_.v;
        let d0 = s0_.v.dot(&e);
        let d1 = e.sqr_magnitude();
        self.closest_a = Float3::lerp_clamp(&polya_[s0_.a], &polya_[s1_.a], d0, d1);
        self.closest_b = Float3::lerp_clamp(&polyb_[s0_.b], &polyb_[s1_.b], d0, d1);
        self.mtv = self.closest_a - self.closest_b;
        self.normalize();
    }

    fn mtv_from_face_case(&mut self, n_:&Float3, polya_:&[Float3], polyb_:&[Float3]) {
        self.mtv = *n_;
        let mut e0 = self.s1.v - self.s0.v;
        let n0 = Float3::cross(&n_, &(*n_-self.s2.v));
        let d0 = n0.dot(&(self.s0.v-self.s2.v));
        let d1 = -n0.dot(&e0);
        let e1 = Float3::lerp(&polya_[self.s0.a], &polya_[self.s1.a], d0, d1);
        let e2 = Float3::lerp(&polyb_[self.s0.b], &polyb_[self.s1.b], d0, d1);
        e0.scale(d0, d1);
        let p = self.s0.v + e0;
        let n0 = p - self.s2.v;
        let d0 = p.dot(&n0);
        let d1 = n0.sqr_magnitude();
        self.closest_a = Float3::lerp(&e1, &polya_[self.s2.a], d0, d1);
        self.closest_b = Float3::lerp(&e2, &polya_[self.s2.b], d0, d1);
        self.normalize();
    }
}