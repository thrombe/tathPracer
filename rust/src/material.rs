
use rand::distributions::{Uniform, Distribution};

use super::vec3d::Vec3d;

use super::ray::Ray;
use super::world::Rng;
use super::math;

#[derive(Debug, Clone)]
pub enum Material {
    Lambertian(Lambertian),
    Metal(Metal),
    Dielectric(Dielectric),
    Lit(Lit), // this can be used as lit if color is greater than 1 and as flat colored stuff if color < 1 (black body if color = 0)
}

#[derive(Debug, Clone)]
pub struct Lambertian {
    pub color: Vec3d,
}

#[derive(Debug, Clone)]
pub struct Metal {
    pub color: Vec3d,
    pub fuzz: f64,
}

#[derive(Debug, Clone)]
pub struct Dielectric {
    pub color: Vec3d,
    pub refractive_index: f64,
    pub fuzz: f64,
}

#[derive(Debug, Clone)]
pub struct Lit {
    pub color: Vec3d,
}

impl Material {
    pub fn scatter(&self, ray: &Ray, normal: &Vec3d, rng: &mut Rng) -> Option<Ray> {
        let rand = Uniform::new(-1.0, 1.0);
        match self {
            Material::Lambertian(_) => {
                let dir = *normal + Vec3d::new(rand.sample(rng), rand.sample(rng), rand.sample(rng));
                // let tea = 0.00000001;
                // if dir.x < tea && dir.y < tea && dir.z < tea { // ray dir goes 0
                //     dir = *normal;
                // }

                Some(Ray::new(ray.pos, dir))
            },
            Material::Metal(obj) => {
                let a = ray.dir.unit();
                let mut reflected = a - *normal*2.0*a.dot(*normal);
                if reflected.dot(*normal) < 0.0 {return None}
                reflected += Vec3d::new(rand.sample(rng), rand.sample(rng), rand.sample(rng))*obj.fuzz;
                Some(Ray::new(ray.pos, reflected))
            },
            Material::Dielectric(obj) => {
                let mut normal = normal.clone();
                let ray_dir_unit = ray.dir.unit();
                let mut ray_dir_dot_normal = ray_dir_unit.dot(normal);

                let ri_ratio = { // ri_matrial2/ri_material1
                    // ray dot normal < 0 means that ray is hitting on the outer surface
                    if ray_dir_dot_normal < 0.0 {1.0/obj.refractive_index}
                    else {
                        ray_dir_dot_normal *= -1.0;
                        normal *= -1.0;
                        obj.refractive_index
                    }

                };
                let cos_i = {
                    if -ray_dir_dot_normal < 1.0 {-ray_dir_dot_normal} else {1.0}
                };
                let sin_sq_i = 1.0-cos_i*cos_i;

                let fuzz = Vec3d::new(rand.sample(rng), rand.sample(rng), rand.sample(rng))*obj.fuzz;

                let schlick_approximation = { // something related to reflection off of glass surfaces at extreme angles
                    let mut r0 = (1.0-ri_ratio)/(1.0+ri_ratio);
                    r0 = r0*r0;
                    let reflectance = r0 + (1.0-r0)*((1.0-cos_i).powf(5.0));
                    reflectance > (rand.sample(rng)+1.0)/1.0
                };
                if ri_ratio * ri_ratio * sin_sq_i > 1.0 || schlick_approximation { // total internal reflection
                
                // if ri_ratio * ri_ratio * sin_sq_i > 1.0 { // total internal reflection
                        Some(Ray::new(ray.pos, ray_dir_unit - normal*2.0*ray_dir_unit.dot(normal) + fuzz))
                } else { // refraction
                    let ray_dir_perp = (ray_dir_unit + normal * cos_i) * ri_ratio;
                    let ray_dir_parallel = normal * (-1.0) * math::abs(1.0 - ray_dir_perp.dot(ray_dir_perp)).sqrt();
                    Some(Ray::new(ray.pos, ray_dir_perp + ray_dir_parallel + fuzz))
                }
            },
            Material::Lit(_) => {
                None
            },
        }
    }

    #[inline(always)]
    pub fn color(&self) -> &Vec3d {
        match self {
            Material::Lambertian(obj) => &obj.color,
            Material::Metal(obj) => &obj.color,
            Material::Dielectric(obj) => &obj.color,
            Material::Lit(obj) => &obj.color,
        }
    }
}

