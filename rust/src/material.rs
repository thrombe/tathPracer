
use rand::distributions::{Uniform, Distribution};

use super::vec3d::Vec3d;

use super::ray::Ray;
use super::world::Rng;


pub enum Material {
    Metal,
    Dielectric,
    Lambertian,
}

impl Material {
    pub fn scatter(&self, ray: &Ray, normal: &Vec3d, rng: &mut Rng) -> Option<Ray> {
        let rand = Uniform::new(-1.0, 1.0);
        match self {
            Material::Lambertian => {
                let dir = *normal + Vec3d::new(rand.sample(rng), rand.sample(rng), rand.sample(rng));
                // let tea = 0.00000001;
                // if dir.x < tea && dir.y < tea && dir.z < tea { // ray dir goes 0
                //     dir = *normal;
                // }

                Some(Ray::new(ray.pos, dir))
            },
            Material::Metal => {
                let a = ray.dir.unit();
                let reflected = a - *normal*2.0*a.dot(*normal);
                if normal.dot(reflected) < 0.0 {return None}
                Some(Ray::new(ray.pos, reflected))
            },
            Material::Dielectric => {
                None
            },
        }
    }
}

