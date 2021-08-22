
use rand::distributions::{Uniform, Distribution};

use super::vec3d::Vec3d;

use super::sphere::Sphere;
use super::material::Material;


pub fn gen_objects() -> Vec<Sphere> {
    let mut objects = Vec::<Sphere>::new();
    objects.push( // Ground
        Sphere {
            center: Vec3d::new(0.0, -1000.0, 0.0),
            radius: 1000.0,
            material: Material::Lambertian,
            color: Vec3d::new(0.5, 0.5, 0.5),
        }
    );

    objects.push(
        Sphere {
            center: Vec3d::new(0.0, 1.0, -10.0),
            radius: 1.0,
            material: Material::Metal,
            color: Vec3d::new(0.77, 1.0, 0.77),
        }
    );
    objects.push(
        Sphere {
            center: Vec3d::new(-2.0, 1.0, -10.0),
            radius: 1.0,
            material: Material::Lambertian,
            color: Vec3d::new(0.65, 1.0, 0.32),
        }
    );
    objects.push(
        Sphere {
            center: Vec3d::new(2.0, 1.0, -10.0),
            radius: 1.0,
            material: Material::Lambertian,
            color: Vec3d::new(0.44, 0.21, 1.0),
        }
    );

    let mut rng = rand::thread_rng();
    let random = Uniform::new(-1.0, 1.0);
    for _ in 0..70 {
        let z = -(random.sample(&mut rng)+1.0)*10.0;
        let x = (random.sample(&mut rng))*10.0;
        objects.push(
            Sphere {
                center: (Vec3d::new(x, 0.0, z) - objects[0].center).unit()*(objects[0].radius+0.4) + objects[0].center,
                radius: 0.4,
                material: Material::Lambertian,
                color: Vec3d::new(random.sample(&mut rng), random.sample(&mut rng), random.sample(&mut rng)),
            }
        );
    }

    objects
}
