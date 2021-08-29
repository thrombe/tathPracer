

use super::vec3d::Vec3d;

use super::material::Material;
use super::world::Rng;

#[derive(Clone)]
pub struct Ray {
    pub pos: Vec3d,
    pub dir: Vec3d,
}

impl Ray {
    pub fn new(pos: Vec3d, dir: Vec3d) -> Self {
        Ray {pos, dir}
    }

    #[inline(always)]
    pub fn new_pos(&mut self, t: f64) {
        self.pos = self.pos + self.dir*t;
    }
}

pub struct RayHitfo { // hit-info?
    pub t: f64,
    pub normal: Vec3d,
    pub material: Material,
    pub ray: Ray,
}

impl RayHitfo {
    #[inline(always)]
    pub fn scatter(&self, rng: &mut Rng) -> Option<Ray> {
        self.material.scatter(&self.ray, &self.normal, rng)
    }
}
