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
mod voxel_octree;
mod object_octree;
mod triangle_octree;
mod aabb;
mod obj_importing;

fn main() {
    let now = time::Instant::now();

    world::run_world();
    // test();
    // scene::gen_world();
    
    println!("{:?}", now.elapsed());
}

pub fn test() {
    // use std::mem::size_of;
    // println!("{:?}", size_of::<Vec<f64>>());
    // use vec3d::Vec3d;
    // use triangle_octree::TriangleOctree;
    // let mut oct = TriangleOctree::new(Vec3d::zero(), 100.0);
    // oct.import_from_obj("../../0builds/objects/axis.obj", 0);
    // dbg!(oct);
}
