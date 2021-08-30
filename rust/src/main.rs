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

    // world::run_world();
    test();
    
    println!("{:?}", now.elapsed());
}

pub fn test() {
    use vec3d::Vec3d;
    use ray::Ray;
    use octree::{Octree};
    use std::mem::size_of;
    // println!("{:?}", size_of::<Vec<f64>>());
    let mut oct = Octree::new(2.0);
    let point = Vec3d::new(0.75, -0.75, -0.75);
    oct.insert_voxel(point, 3);
    // println!("{:?}", oct);
    println!("");
    let ray = Ray::new(Vec3d::new(-1.001, 1.002, 1.003), Vec3d::new(1.0, -1.0, -1.0));
    println!("{:?}", ray);
    let hitfo = oct.hit(&ray, 0.00001);
    println!("{:?}", hitfo);
    
}
