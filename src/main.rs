mod vector;
mod minkowski_sum_helper;
mod gjk3d;
use vector::Float3;
use std::fs::File;
use std::io::{self, BufRead};

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
    match read_input("input") {
        Ok((polya, polyb)) => match gjk3d::check(&polya, &polyb, &Float3::zero()) {
            Some((collide, info)) => { println!("collide"); }
            None => { println!("error"); }
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

fn read_input(path:&str) -> Result<(Vec<Float3>, Vec<Float3>), ParseError> {
    let file = File::open(path).map_err(ParseError::IO)?;
    let mut iter = io::BufReader::new(file).lines().filter_map(|line| line.ok()).map(|line| parse_from(&line));
    let a = iter.next().ok_or(ParseError::Length)??;
    let b = iter.next().ok_or(ParseError::Length)??;
    Ok((a, b))
}

fn parse_from(line: &str) -> Result<Vec<Float3>, ParseError> {
    let mut iter = line.split(|c| c == ',' || c == ' ').filter_map(|s| s.parse::<f32>().ok());
    let mut ret:Vec<Float3> = Vec::new();
    let mut x = iter.next().ok_or(ParseError::Non3Multiple)?;
    loop {
        ret.push(Float3 {
            x:x,
            y:iter.next().ok_or(ParseError::Non3Multiple)?,
            z:iter.next().ok_or(ParseError::Non3Multiple)?,
            });
        match iter.next() {
            Some(v) => x = v,
            None => break,
        }
    }
    Ok(ret)
}