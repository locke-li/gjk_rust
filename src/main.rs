mod vector;
mod simplex_based_cd_helper;
mod gjk3d;
use vector::Float3;
use std::fs::File;
use std::io::{self, BufRead};
type Polygon = (Vec<Float3>, Float3);

fn main() {
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
    let mut frame = simplex_based_cd_helper::Frame3Simplex::<Float3>::new();
    match read_input("input") {
        Ok(((poly_a, center_a), (poly_b, center_b))) => 
            match gjk3d::check(&poly_a, &poly_b, &(center_a - center_b), &mut frame) {
                Ok(true) => { println!("collide {:?}", frame); }
                Ok(false) => { println!("miss {:?}", frame); }
                Err(e) => { println!("error {:?}", e); }
            },
        Err(s) => println!("{:?}", s),
    }
}

#[derive(Debug)]
enum ParseError {
    IO(io::Error),
    Parse(std::num::ParseFloatError),
    Length,
    Non3Multiple,
}

fn read_input(path:&str) -> Result<(Polygon, Polygon), ParseError> {
    let file = File::open(path).map_err(ParseError::IO)?;
    let mut iter = io::BufReader::new(file).lines().filter_map(|line| line.ok()).map(|line| parse_from(&line));
    let a = iter.next().ok_or(ParseError::Length)?.ok_or(ParseError::Non3Multiple)?;
    let b = iter.next().ok_or(ParseError::Length)?.ok_or(ParseError::Non3Multiple)?;
    Ok((a, b))
}

fn parse_from(line: &str) -> Option<Polygon> {
    let mut iter = line.split(|c| c == ',' || c == ' ').filter_map(|s| s.parse::<f32>().ok());
    let mut ret:Vec<Float3> = Vec::new();
    let pos_x = iter.next()?;
    let pos_y = iter.next()?;
    let pos_z = iter.next()?;
    let mut center_x = 0.0;
    let mut center_y = 0.0;
    let mut center_z = 0.0;
    let mut x = iter.next()?;
    loop {
        let y = iter.next()?;
        let z = iter.next()?;
        center_x += x;
        center_y += y;
        center_z += z;
        ret.push(Float3{x:x + pos_x, y:y + pos_y, z:z + pos_z});
        match iter.next() {
            Some(v) => x = v,
            None => break,
        }
    }
    let len = ret.len() as f32;
    Some((ret, Float3{x:center_x/len + pos_x, y:center_y/len + pos_y, z:center_z/len + pos_z}))
}