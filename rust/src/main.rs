#![allow(dead_code)]
#![allow(non_snake_case)]

use std::time;
mod img;
// mod math;
mod vec3d;
mod progress_indicator;

fn main() {
    let now = time::Instant::now();

    
    println!("{:?}", now.elapsed());
}
