#![allow(dead_code)]
#![allow(non_snake_case)]

use std::time;
mod img;
// mod math;
mod vec3d;
mod progress_indicator;
mod world;

fn main() {
    let now = time::Instant::now();

    world::run_world();
    // world::test();
    
    println!("{:?}", now.elapsed());
}
