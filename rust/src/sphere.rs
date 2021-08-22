

use super::vec3d::Vec3d;

use super::material::Material;
use super::ray::Ray;
use super::world::Rng;


pub struct Sphere {
    pub center: Vec3d,
    pub radius: f64,
    pub material: Material,
}

impl Sphere {
    /// returns t for the ray if hit, else returns None
    pub fn hit(&self, ray: &Ray, t_correction: f64) -> Option<f64> {
        let oc = ray.pos-self.center;
        let neg_b = -ray.dir.dot(oc);
        let b_sq = ray.dir.dot(ray.dir);
        let d_by_4 = neg_b*neg_b - b_sq*(oc.dot(oc) - self.radius*self.radius);
        if d_by_4 < 0.0 { // didnt hit
            return None
        }
        let t = (neg_b - d_by_4.sqrt())/b_sq;
        if t < t_correction { // ray hitting behind the camera or really close to object ( t < 0.0000001 for really close thing)
            return None
        }
        Some(t)
    }

    pub fn scatter(&self, ray: &Ray, rng: &mut Rng) -> Option<Ray> {
        self.material.scatter(ray, &mut self.normal(&ray.pos), rng)
    }

    pub fn normal(&self, point: &Vec3d) -> Vec3d {
        (*point - self.center).unit()
    }

    pub fn color(&self) -> &Vec3d {
        self.material.color()
    }
}
