mod vector;
mod minkowski_sum_helper;
mod gjk3d;
use vector::Float3;

fn main() {
    let a = Float3 {x:1.0, y:1.0, z:1.0};
    let b = Float3 {x:2.0, y:3.0, z:4.0};
    let a_copy = Float3 {x:1.0, y:1.0, z:1.0};
    println!("{} {}", a == b, a == a_copy);
    println!("{}", a.sqr_magnitude());
    let c = a + b;
    println!("{:?}", c);
    let r = a.dot(&b);
    println!("a dot b = {}", r);
    println!("a cross b = {:?}", Float3::cross(&a, &b));
}
