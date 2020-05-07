mod vector;
mod simplex_based_cd_helper;
mod gjk3d;
use vector::Float3;
use std::fs::File;
use std::io::{self, BufRead};

struct ObjectInput {
    pub vertex: Vec<Float3>,
    pub position: Float3,
    pub center: Float3,
    pub velocity: Float3,
}

fn main() {
    vector_test();
    let mut frame = simplex_based_cd_helper::Frame3Simplex::<Float3>::new();
    //TODO add main loop to animate the objects
    match read_input("input") {
        Ok((obj_a, obj_b)) => 
            match gjk3d::check(&obj_a.vertex, &obj_b.vertex, &(obj_a.center - obj_b.center), &mut frame) {
                Ok(true) => { println!("collide\n{:?}", frame); }
                Ok(false) => { println!("miss\n{:?}", frame); }
                Err(e) => { println!("error\n{:?}", e); }
            },
        Err(s) => println!("{:?}", s),
    }
}

#[derive(Debug)]
enum ParseError {
    IO(io::Error),
    Length,
    Non3Multiple,
}

fn vector_test() {
    let a = Float3 {x:1.0, y:1.0, z:1.0};
    let b = Float3 {x:2.0, y:3.0, z:4.0};
    let a_copy = Float3 {x:1.0, y:1.0, z:1.0};
    println!("{} {}", a == b, a == a_copy);
    println!("{}", a.sqr_magnitude());
    let c = a + b;
    println!("{}", c);
    let r = a.dot(&b);
    println!("a dot b = {}", r);
    println!("a cross b = {}", Float3::cross(&a, &b));
}

fn read_input(path:&str) -> Result<(ObjectInput, ObjectInput), ParseError> {
    let file = File::open(path).map_err(ParseError::IO)?;
    let mut iter = io::BufReader::new(file).lines().filter_map(|line| line.ok())
        .filter_map(|line| if line.starts_with("#") { Option::None } else { Option::Some(line) })
        .map(|line| parse_from(&line));
    let a = iter.next().ok_or(ParseError::Length)?.ok_or(ParseError::Non3Multiple)?;
    let b = iter.next().ok_or(ParseError::Length)?.ok_or(ParseError::Non3Multiple)?;
    Ok((a, b))
}

fn parse_from(line:&str) -> Option<ObjectInput> {
    let mut iter = line.split(|c| c == ',' || c == '|' || c == ' ').filter_map(|s| s.parse::<f32>().ok());
    let mut ret:Vec<Float3> = Vec::new();
    //root/origin position
    let pos_x = iter.next()?;
    let pos_y = iter.next()?;
    let pos_z = iter.next()?;
    //move velocity
    let vel_x = iter.next()?;
    let vel_y = iter.next()?;
    let vel_z = iter.next()?;
    let mut center_x = 0.0;
    let mut center_y = 0.0;
    let mut center_z = 0.0;
    //x component of the first input vertex, at least one point is required
    let mut x = iter.next()?;
    loop {
        let y = iter.next()?;
        let z = iter.next()?;
        center_x += x;
        center_y += y;
        center_z += z;
        ret.push(Float3{x:x + pos_x, y:y + pos_y, z:z + pos_z});
        //next x
        match iter.next() {
            Some(v) => x = v,
            None => break,
        }
    }
    let len = ret.len() as f32;
    Some(ObjectInput {
        vertex: ret,
        position: Float3{x:pos_x, y:pos_y, z:pos_z},
        center: Float3{x:center_x/len + pos_x, y:center_y/len + pos_y, z:center_z/len + pos_z}, 
        velocity: Float3{x:vel_x, y:vel_y, z:vel_z}
    })
}