

use super::vec3d::Vec3d;

pub struct Aabb { // axis aligned bounding box
    pub min: Vec3d,
    pub max: Vec3d,
}

impl Aabb {
    pub fn new(min: Vec3d, max: Vec3d) -> Self {
        Self {min, max}
    }

    pub fn new_from_points(points: Vec<Vec3d>) -> Self {
        let mut min = Vec3d::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
        let mut max = Vec3d::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
        for point in points {
            if point.x > max.x {max.x = point.x}
            if point.y > max.y {max.y = point.y}
            if point.z > max.z {max.z = point.z}
            if point.x < min.x {min.x = point.x}
            if point.y < min.y {min.y = point.y}
            if point.z < min.z {min.z = point.z}
        }
        Self::new(min, max)
    }

    pub fn contains(&self, other: &Self) -> bool {
        if (other.min.x < self.min.x) || (other.min.y < self.min.y) || (other.min.z < self.min.z) || 
           (other.max.x > self.max.x) || (other.max.y > self.max.y) || (other.max.z > self.max.z) {return false}
        true
    }
}

use std::ops::{Add, Sub, Mul};
impl Add<Vec3d> for Aabb {
    type Output = Self;
    #[inline(always)]
    fn add(mut self, vec: Vec3d) -> Self{
        self.min += vec;
        self.max += vec;
        self
    }
}
impl Sub<Vec3d> for Aabb {
    type Output = Self;
    #[inline(always)]
    fn sub(mut self, vec: Vec3d) -> Self {
        self.min -= vec;
        self.max -= vec;
        self
    }
}
impl Mul<f64> for Aabb {
    type Output = Self;
    #[inline(always)]
    fn mul(mut self, t: f64) -> Self {
        self.min *= t;
        self.max *= t;
        self
    }
}
