

use super::vec3d::Vec3d;
// use super::math;

use super::material::Material;
use super::ray::Ray;
use super::world::Rng;

pub enum Object {
    Sphere(Sphere),
    Plane(Plane),
    Triangle(Triangle),
}

impl Object {
    pub fn hit(&self, ray: &Ray, t_correction: f64) -> Option<f64> {
        match self {
            Object::Sphere(obj) => obj.hit(ray, t_correction),
            Object::Plane(obj) => obj.hit(ray, t_correction),
            Object::Triangle(obj) => obj.hit(ray, t_correction),
        }
    }

    pub fn scatter(&self, ray: &Ray, rng: &mut Rng) -> Option<Ray> {
        match self {
            Object::Sphere(obj) => obj.scatter(ray, rng),
            Object::Plane(obj) => obj.scatter(ray, rng),
            Object::Triangle(obj) => obj.scatter(ray, rng),
        }    
    }

    pub fn normal(&self, ray: &Ray) -> Vec3d {
        match self {
            Object::Sphere(obj) => obj.normal(ray),
            Object::Plane(obj) => obj.normal(ray),
            Object::Triangle(obj) => obj.normal(),
        }
    }

    pub fn color(&self) -> &Vec3d {
        match self {
            Object::Sphere(obj) => obj.color(),
            Object::Plane(obj) => obj.color(),
            Object::Triangle(obj) => obj.color(),
        }
    }
    
    pub fn material(&self) -> &Material {
        match self {
            Object::Sphere(obj) => &obj.material,
            Object::Plane(obj) => &obj.material,
            Object::Triangle(obj) => &obj.material,
        }
    }
}

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
        let t = { // minimum t when ray originates outside, min +ve t when ray originates inside
            let sqrt_d_by_4 = d_by_4.sqrt();
            let m_t = (neg_b - sqrt_d_by_4)/b_sq; // when outside sphere, m_t only is enough
            let p_t = (neg_b + sqrt_d_by_4)/b_sq; // when inside sphere, p_t is good (i think)

            // if p_t -ve, {m_t is +ve -> good, m_t -ve -> then the t_correction cancels it}
            // if m_t -ve, we know p_t is +ve -> good
            // if both +ve -> min(both vals)
            if p_t < 0.0 {m_t} else if m_t < 0.0 {p_t} else if p_t < m_t {p_t} else {m_t}
        };
        if t < t_correction { // ray hitting behind the camera or really close to object ( t < 0.0000001 for really close thing)
            return None
        }
        Some(t)
    }

    pub fn scatter(&self, ray: &Ray, rng: &mut Rng) -> Option<Ray> {
        self.material.scatter(ray, &mut self.normal(ray), rng)
    }

    pub fn normal(&self, ray: &Ray) -> Vec3d {
        (ray.pos - self.center).unit()
    }

    pub fn color(&self) -> &Vec3d {
        self.material.color()
    }
}

pub struct Plane {
    pub normal: Vec3d,
    pub point: Vec3d,
    pub material: Material,
}

impl Plane {
    pub fn hit(&self, ray: &Ray, t_correction: f64) -> Option<f64> {
        // (a+bt-p).n = 0 -> t = ((p-a).n)/(b.n) : ray is a+bt, p is a point on plane, n is normal to plane

        // let normal = self.normal(ray);
        let normal = self.normal;
        let t = (self.point-ray.pos).dot(normal)/ray.dir.dot(normal); // normal or -ve normal dosent matter here cuz normal is in both num and denom
        if t > t_correction {Some(t)} else {None}
    }

    pub fn scatter(&self, ray: &Ray, rng: &mut Rng) -> Option<Ray> {
        self.material.scatter(ray, &mut self.normal(ray), rng)
    }

    pub fn normal(&self, ray: &Ray) -> Vec3d {
        if self.normal.dot(ray.dir) < 0.0 {self.normal.unit()} else {self.normal.unit()*(-1.0)}
    }

    pub fn color(&self) -> &Vec3d {
        self.material.color()
    }
}

pub struct Triangle {
    pub vertices: (Vec3d, Vec3d, Vec3d),
    pub normal: Vec3d,
    pub material: Material,
}

impl Triangle {
    pub fn hit(&self, ray: &Ray, t_correction: f64) -> Option<f64> {
        // M-1
        // let OAB (clockwise) be a triangle if A.cross(P).dot(N) > 0 and B.cross(P).dot(N) < 0 -> P is inside
        // N is normal = A.cross(B), and P is the required point

        // M-2
        // let ABC be triangle, a+bt be ray, a+bt=A(1-w2-w3)+Bw2+Cw3 (barycentric coords)
        // A-a=[b, A-B, A-C]*[[t], [w2], [w3]] then apply crammers rule
        // (remember b, A-B .. are vectors, so that is a 3x3 matrix)
        // calculate D, D1, D2, D3, if D != 0(close to 0 cuz float)(if D<0, triangle is backwards), t = D1/D ....
        // D can be calculated by a.dot(b.cross(c))

        // finding t by plane intersection
        let t = (self.vertices.0-ray.pos).dot(self.normal)/ray.dir.dot(self.normal);
        if t < t_correction {return None}

        let A = self.vertices.0;
        let B = self.vertices.1;
        let C = self.vertices.2;
        let a = ray.pos;
        let b = ray.dir;
        let Ama = A-a;
        let AmB = A-B;
        let AmC = A-C;
        let D = b.dot(AmB.cross(AmC));
        if !(D > t_correction || D < -t_correction) {return None}
        // let t = Ama.dot(AmB.cross(AmC))/D; // D1/D
        // if t < t_correction {return None}
        let w2 = b.dot(Ama.cross(AmC))/D; // D2/D
        let w3 = b.dot(AmB.cross(Ama))/D; // D3/D
        if !(w2 > 0.0 && w3 > 0.0 && w2+w3 < 1.0) {return None}
        Some(t)
    }

    pub fn scatter(&self, ray: &Ray, rng: &mut Rng) -> Option<Ray> {
        self.material.scatter(ray, &mut self.normal(), rng)
    }

    pub fn normal(&self) -> Vec3d {
        self.normal // it has a defined outer direction (vertices clockwise/anti-clockwise)
    }

    pub fn color(&self) -> &Vec3d {
        self.material.color()
    }
}
