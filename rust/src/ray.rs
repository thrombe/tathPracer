


use super::vec3d::Vec3d;


pub struct Ray {
    pub pos: Vec3d,
    pub dir: Vec3d,
}

impl Ray {
    pub fn new(pos: Vec3d, dir: Vec3d) -> Self {
        Ray {pos, dir}
    }

    pub fn at(&mut self, t: f64) {
        self.pos = self.pos + self.dir*t;
    }
}
