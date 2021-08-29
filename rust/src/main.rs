#![allow(dead_code)]
#![allow(non_snake_case)]

use std::time;

mod img;
mod math;
mod vec3d;
mod progress_indicator;

mod world;
mod ray;
mod material;
mod objects;
mod scene;
mod octree;

fn main() {
    let now = time::Instant::now();

    world::run_world();
    // test();
    
    println!("{:?}", now.elapsed());
}

pub fn test() {
    use vec3d::Vec3d;
    use octree::{Octree};
    use std::mem::size_of;
    println!("{:?}", size_of::<Vec<f64>>());
    let mut oct = Octree::new(8.0);
    let point = Vec3d::new(0.3, 2.45, 1.66);
    oct.insert_voxel(point, 0);
    println!("{:?}", oct);
}
