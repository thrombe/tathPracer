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

fn main() {
    let now = time::Instant::now();

    world::run_world();
    // test();
    
    println!("{:?}", now.elapsed());
}

pub fn test() {
    use objects::Sphere;
    use material::{Material, Lit};
    use vec3d::Vec3d;
    use ray::Ray;
    
    let sphere = Sphere {
            center: Vec3d::new(0.0, 0.0, 0.0),
            radius: 500.0,
            material: Material::Lit(
                Lit {
                    color: Vec3d::new(0.8, 3.7, 0.9),
                }
            ),
        };

    let ray = Ray {
        pos:Vec3d::new(0.0, 1.0, 0.0),
        dir:Vec3d::new(0.0, 0.0, -1.0),
    };
    let hit = sphere.hit(&ray, 0.0001);
    println!("{:?}", hit);
}