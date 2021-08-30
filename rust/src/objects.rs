
use std::sync::Arc;

use super::vec3d::Vec3d;
// use super::math;

use super::material::{Material, Lit};
use super::ray::{Ray, RayHitfo};
use super::world::Rng;
use super::octree::{Octree, OctreeBranch, OctreePos};

pub enum Object {
    Sphere(Sphere),
    Plane(Plane),
}

impl Object {
    pub fn hit(&self, ray: &Ray, t_correction: f64) -> Option<RayHitfo> {
        match self {
            Object::Sphere(obj) => obj.hit(ray, t_correction),
            Object::Plane(obj) => obj.hit(ray, t_correction),
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
    pub fn hit(&self, ray: &Ray, t_correction: f64) -> Option<RayHitfo> {
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
        
        let mut ray = ray.clone();
        ray.new_pos(t);
        Some(RayHitfo {
            t,
            normal: self.normal(&ray),
            ray,
            material: self.material.clone(),
        })
    }

    #[inline(always)]
    pub fn normal(&self, ray: &Ray) -> Vec3d {
        (ray.pos - self.center).unit()
    }
}

pub struct Plane {
    pub normal: Vec3d,
    pub point: Vec3d,
    pub material: Material,
}

impl Plane {
    pub fn hit(&self, ray: &Ray, t_correction: f64) -> Option<RayHitfo> {
        // (a+bt-p).n = 0 -> t = ((p-a).n)/(b.n) : ray is a+bt, p is a point on plane, n is normal to plane

        let normal = self.normal;
        let t = (self.point-ray.pos).dot(normal)/ray.dir.dot(normal); // normal or -ve normal dosent matter here cuz normal is in both num and denom
        if t > t_correction {
            let mut ray = ray.clone();
            ray.new_pos(t);
            Some(RayHitfo {
                t,
                normal: self.normal(&ray),
                ray,
                material: self.material.clone(),
            })
        } else {None}
    }

    #[inline(always)]
    pub fn normal(&self, ray: &Ray) -> Vec3d {
        if self.normal.dot(ray.dir) < 0.0 {self.normal.unit()} else {self.normal.unit()*(-1.0)}
    }
}

impl Octree {
    pub fn hit(&self, ray: &Ray, t_correction: f64) -> Option<RayHitfo> {
        // intersect the bb and return if not hit
        // todo

        let ray = Ray::new(self.world_to_tree_space(ray.pos), ray.dir);
        
        // 1
        let t = |P: Vec3d| { // sorted by min
            let (mut r, mut u, mut b) = OctreePos::get_rub_masks();
            let mut ray = ray.clone();
            // flipping the ray wrt the mid planes (local)
            // note that the ray isnt actually being used in the rest of the hit code. so no need to modify it in other places
            if ray.dir.x < 0.0 {
                ray.dir.x *= -1.0;
                ray.pos.x *= -1.0;
                // ray.pos.x = 2.0-ray.pos.x;
                r = !r;
            }
            if ray.dir.y < 0.0 {
                ray.dir.y *= -1.0;
                ray.pos.y *= -1.0;
                // ray.pos.y = 2.0-ray.pos.y;
                u = !u;
            }
            if ray.dir.z < 0.0 {
                ray.dir.z *= -1.0;
                ray.pos.z *= -1.0;
                // ray.pos.z = 2.0-ray.pos.z;
                b = !b;
            }
            let tx = if ray.dir.x != 0.0 {(P.x-ray.pos.x)/ray.dir.x} else {1000000000000000.0}; // yz - r
            let ty = if ray.dir.y != 0.0 {(P.y-ray.pos.y)/ray.dir.y} else {1000000000000000.0}; // xz - u
            let tz = if ray.dir.z != 0.0 {(P.z-ray.pos.z)/ray.dir.z} else {1000000000000000.0}; // xy - b
            let mut ts = ((tx, r), (ty, u), (tz, b));
            if ts.0.0 > ts.1.0 {
                let temp = ts.0;
                ts.0 = ts.1;
                ts.1 = temp;
            }
            if ts.1.0 > ts.2.0 {
                let temp = ts.1;
                ts.1 = ts.2;
                ts.2 = temp;
            }
            if ts.0.0 > ts.1.0 {
                let temp = ts.0;
                ts.0 = ts.1;
                ts.1 = temp;
            }
            ts
        };
        let t0 = t(Vec3d::new(-1.0, -1.0, -1.0));
        let t1 = t(Vec3d::new(1.0, 1.0, 1.0));
        let (tm, del_t) = {
            let mut tm = t0;
            let del_t = ((t1.0.0-t0.0.0)/2.0, (t1.1.0-t0.1.0)/2.0, (t1.2.0-t0.2.0)/2.0);
            tm.0.0 += del_t.0;
            tm.1.0 += del_t.1;
            tm.2.0 += del_t.2;
            (tm, del_t)
        };



        if let Some(hitfo) = self.main_branch.hit(&ray, t0.2, del_t, tm) {
           Some(RayHitfo {
            //    t: self.tree_to_world_space_f64(hitfo.t),
               ray: Ray::new(self.tree_to_world_space(hitfo.ray.pos), hitfo.ray.dir),
               ..hitfo
           })
        } else {
            None
        }
    }
}

impl OctreeBranch {
    pub fn hit(&self, ray: &Ray, t0: (f64, u8), del_t: (f64, f64, f64), tm: ((f64, u8), (f64, u8), (f64, u8))) 
        -> Option<RayHitfo> {
        let mut ray2 = ray.clone();//////
        ray2.new_pos(t0.0);////////
        println!("{:?}, t {:?}, raypos {:?}, tm {:?}", self.pos, t0, ray2.pos, tm);//////
        println!("{:?} t1-> {:?}", t0, 0);
        let masks = (
            (t0.0, !t0.1), (tm.0.0, if tm.0.0 < t0.0 {tm.0.1} else {!tm.0.1}),
            (tm.1.0, if tm.1.0 < t0.0 {tm.1.1} else {!tm.1.1}),
            (tm.2.0, if tm.2.0 < t0.0 {tm.2.1} else {!tm.2.1})
        );
        let entry_child_mask = masks.0.1 & masks.1.1 & masks.2.1;
        println!("passing through ecm {:?}, t {:?}", OctreePos::new(entry_child_mask), t0.0);
        if entry_child_mask & self.branch_mask > 0 {
            // visit the child
            // println!("ecm {:?}", entry_child_mask);////////
            let child = self.get_branch(OctreePos::new(entry_child_mask));
            let del_t = (del_t.0/2.0, del_t.1/2.0, del_t.2/2.0);
            let mut tm = tm;
            tm.0.0 -= del_t.0;
            tm.1.0 -= del_t.1;
            tm.2.0 -= del_t.2;
            if let Some(hitfo) = child.hit(ray, t0, del_t, tm) {
                return Some(hitfo)
            }
        } else if entry_child_mask & self.child_mask > 0 {
            // its a voxel. so return?
            let t = t0.0;
            let mut ray = ray.clone();
            ray.new_pos(t);
            return Some(RayHitfo {
                t,
                normal: Vec3d::new(1.0, 1.0, 1.0),
                material: Material::Lit(Lit {color: Vec3d::zero()}),
                ray,
            })
        }
        // its an off voxel or its inner voxels were not hit. so visit next
        let mouve = |child: u8, plane_mask: u8| -> Option<u8> {
            if child & plane_mask > 0 {return None} // child is already in that side (eg -> move right (but im already in right, ig ill exit))
            let mut other_masks = {
                // if plane_mask == masks.1.1 || plane_mask == !masks.1.1 {(masks.2.1, masks.3.1)}
                // else if plane_mask == masks.2.1 || plane_mask == !masks.2.1 {(masks.1.1, masks.3.1)}
                // else {(masks.1.1, masks.2.1)}
                if plane_mask == tm.0.1 || plane_mask == !tm.0.1 {(tm.1.1, tm.2.1)}
                else if plane_mask == tm.1.1 || plane_mask == !tm.1.1 {(tm.0.1, tm.2.1)}
                else {(tm.0.1, tm.1.1)}
            };
            other_masks.0 = if other_masks.0 & child > 0 {other_masks.0} else {!other_masks.0};
            other_masks.1 = if other_masks.1 & child > 0 {other_masks.1} else {!other_masks.1};
            let next = plane_mask & other_masks.0 & other_masks.1;
            println!("move {:?}, plane_mask {:#010b}, other_masks {:?}", next, plane_mask, other_masks);////////
            // println!("move pos {:?}", OctreePos::new(next));////////
            // if next == child {return None}
            Some(next)
        };
        // if tm.0.0 > t1.0.0 {return None}
        let c1: u8;
        match mouve(entry_child_mask, tm.0.1) {
            Some(ueight) => c1 = ueight,
            None => return None,
        }
        println!("passing through c1 {:?}, t {:?}", OctreePos::new(c1), tm.0.0);
        //enter c1
        if c1 & self.branch_mask > 0 {
            // visit the child
            println!("c1 {:?}", c1);////////
            let child = self.get_branch(OctreePos::new(c1));
            let del_t = (del_t.0/2.0, del_t.1/2.0, del_t.2/2.0);
            let t0 = tm.0;
            let mut tm = tm;
            tm.0.0 += del_t.0;
            tm.1.0 -= del_t.1;
            tm.2.0 -= del_t.2;
            if let Some(hitfo) = child.hit(ray, t0, del_t, tm) {
                return Some(hitfo)
            }
        } else if c1 & self.child_mask > 0 {
            // its a voxel. so return?
            let t = tm.0.0;
            let mut ray = ray.clone();
            ray.new_pos(t);
            return Some(RayHitfo {
                t,
                normal: Vec3d::new(1.0, 1.0, 1.0),
                material: Material::Lit(Lit {color: Vec3d::zero()}),
                ray,
            })
        }
        // if tm.1.0 > t1.1.0 {return None}
        let c2: u8;
        match mouve(c1, tm.1.1) {
            Some(ueight) => c2 = ueight,
            None => return None,
        }
        println!("passing through c2 {:?}, t {:?}", OctreePos::new(c2), tm.1.0);
        //enter c2
        if c2 & self.branch_mask > 0 {
            // visit the child
            println!("c2 {:?}, c1 {:?}, tm.1.1 {:?}, masks {:?}", c2, c1, tm.1.1, masks);////////
            let child = self.get_branch(OctreePos::new(c2));
            let del_t = (del_t.0/2.0, del_t.1/2.0, del_t.2/2.0);
            let t0 = tm.1;
            let mut tm = tm;
            tm.0.0 += del_t.0;
            tm.1.0 += del_t.1;
            tm.2.0 -= del_t.2;
            if let Some(hitfo) = child.hit(ray, t0, del_t, tm) {
                return Some(hitfo)
            }
        } else if c2 & self.child_mask > 0 {
            // its a voxel. so return?
            let t = tm.1.0;
            let mut ray = ray.clone();
            ray.new_pos(t);
            return Some(RayHitfo {
                t,
                normal: Vec3d::new(1.0, 1.0, 1.0),
                material: Material::Lit(Lit {color: Vec3d::zero()}),
                ray,
            })
        }
        // if tm.2.0 > t1.2.0 {return None}
        let c3: u8;
        match mouve(c2, tm.2.1) {
            Some(ueight) => c3 = ueight,
            None => return None,
        }
        println!("passing through c3 {:?}, t {:?}", OctreePos::new(c3), tm.2.0);
        //enter c3
        if c3 & self.branch_mask > 0 {
            // visit the child
            // let mut ray2 = ray.clone();
            // ray2.new_pos(tm.2);
            let child = self.get_branch(OctreePos::new(c3));
            let del_t = (del_t.0/2.0, del_t.1/2.0, del_t.2/2.0);
            let t0 = tm.2;
            let mut tm = tm;
            tm.0.0 += del_t.0;
            tm.1.0 += del_t.1;
            tm.2.0 += del_t.2;
            println!("c3 {:?}, tm {:?}", c3, tm);////////
            if let Some(hitfo) = child.hit(ray, t0, del_t, tm) {
                return Some(hitfo)
            }
        } else if c3 & self.child_mask > 0 {
            // its a voxel. so return?
            println!("{:#010b}", self.child_mask);
            let t = tm.2.0;
            let mut ray = ray.clone();
            ray.new_pos(t);
            return Some(RayHitfo {
                t,
                normal: Vec3d::new(1.0, 1.0, 1.0),
                material: Material::Lit(Lit {color: Vec3d::zero()}),
                ray,
            })
        }
        None
    }
}

// refactoring needed for rayhitfo (down)

pub struct Triangle {
    vertex_indices: (u32, u32, u32),
    pub normal: Vec3d,
    pub material: Material,
}

impl Triangle {
    pub fn hit(&self, vertices: (Vec3d, Vec3d, Vec3d), ray: &Ray, t_correction: f64) -> Option<f64> {
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
        let t = (vertices.0-ray.pos).dot(self.normal)/ray.dir.dot(self.normal);
        if t < t_correction {return None}

        let A = vertices.0;
        let B = vertices.1;
        let C = vertices.2;
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

    fn get_vertices(&self, vertices: Arc<Vec<Vec3d>>) -> (Vec3d, Vec3d, Vec3d) {
        (vertices[self.vertex_indices.0 as usize], vertices[self.vertex_indices.1 as usize], vertices[self.vertex_indices.2 as usize])
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

pub struct TriangleMesh {
    vertices: Arc/* or Box?(why rc for something inside a struct)*/<Vec<Vec3d>>, // meshes inside meshes??
    triangles: Arc<Vec<Triangle>>,
    // bounding_box: (Vec3d, Vec3d), // min, max // axis aligned
    position: Vec3d, // we can use this to move objects without editing all coords (notes)
    bounding_box: (Vec3d, Vec3d, Vec3d, Vec3d), // position of a corner, xdir, ydir, zdir // not axis aligned // we can get axis aligned using this
}
// read notes

