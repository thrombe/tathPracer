
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
        
        // does this work when ray originates from side the tree?

        let ray = Ray::new(self.world_to_tree_space(ray.pos), ray.dir);
        let mut planes = OctreePos::get_rub_masks();
        let mut bbox_min = Vec3d::new(-1.0, -1.0, -1.0);
        let mut dt = {
            let mut dt = Vec3d::new(1.0/ray.dir.x, 1.0/ray.dir.y, 1.0/ray.dir.z); // f64::Infinity comes up if any is 0.0
            // either flip the ray wrt the mid planes (local)
            // or choose the correct bbox edges
            if ray.dir.x < 0.0 {
                dt.x *= -1.0;
                bbox_min.x *= -1.0;
                planes.0 = !planes.0;
            }
            if ray.dir.y < 0.0 {
                dt.y *= -1.0;
                bbox_min.y *= -1.0;
                planes.1 = !planes.1;
            }
            if ray.dir.z < 0.0 {
                dt.z *= -1.0;
                bbox_min.z *= -1.0;
                planes.2 = !planes.2;
            }
            dt
        };
        let mut t0 = Vec3d::new((bbox_min.x-ray.pos.x)/dt.x, (bbox_min.y-ray.pos.y)/dt.y, (bbox_min.z-ray.pos.z)/dt.z);
        
        fn swap<T: Copy>(a: &mut T, b: &mut T) {
            let temp = *a;
            *a = *b;
            *b = temp;
        }
        { // sort in ascending order
            if t0.x < t0.y {
                swap(&mut t0.x, &mut t0.y);
                swap(&mut dt.x, &mut dt.y);
                swap(&mut planes.0, &mut planes.1)
            }
            if t0.y < t0.z {
                swap(&mut t0.y, &mut t0.z);
                swap(&mut dt.y, &mut dt.z);
                swap(&mut planes.1, &mut planes.2)
            }
            if t0.x < t0.y {
                swap(&mut t0.x, &mut t0.y);
                swap(&mut dt.x, &mut dt.y);
                swap(&mut planes.0, &mut planes.1)
            }
        }
        let tm = t0 + dt;

        if let Some(hitfo) = self.main_branch.hit(&ray, t0.z, planes.2, planes, &dt, &tm) {
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

impl OctreeBranch { // t, t_plane, planes, dt, tm

    // think about all explanation in this func as if ray.dir.x and y and z > 0
    // the planes are masks of voxels to right/up/down/whatever. so if we "&" multiple of these, we can constrict where the required voxel is
    // eg: let 000 represent x, y, z, so 111 & 010 would mean gimme something (out of 111) that has y = 1 (this makes more sense with u8 but im lazy)
    // note that order of dt, tm and planes is related. i.e planes.0, tm.x, dt.x are for same plane
    pub fn hit(&self, ray: &Ray, t: f64, t_plane: u8, planes: (u8, u8, u8), dt: &Vec3d, tm: &Vec3d) -> Option<RayHitfo> {
        // debugging
        let mut ray2 = ray.clone();//////
        ray2.new_pos(t);////////
        println!("{:?}, t {:?}, raypos {:?}", self.pos, t, ray2.pos);//////

        let dt = *dt*0.5; // dt for every sub voxel is half of its parents (since size is half)

        // to choose what subvoxel to enter, we know, it entered crossing t_plane. so it should be a voxel thats to that side of the main voxel
        // eg: if we cross the left wall(bbox_min), we know the first sub voxel must be touching this wall.
        // since, by default right, up and forward are chosen, left would mean !t_plane
        // similarly, to choose what side of the other mid planes the sub voxel lies on, we can compare if the t value for that plane and the main voxel t value
        let entry_child_mask = !t_plane
                               & if tm.x < t {planes.0} else {!planes.0}
                               & if tm.y < t {planes.1} else {!planes.1};
        // t value for first sub voxel is same as of the main voxel
        if let Some(hitfo) = self.try_hit_subvoxel(entry_child_mask, ray, t, t_plane, planes, &dt, &(*tm-dt)) {
            return Some(hitfo)
        }
        // its an off voxel or its inner voxels were not hit. so visit next

        // crossing the first mid plane (planes.0)
        let c1: u8;
        match self.get_next_voxel(entry_child_mask, planes.0, planes) {
            Some(next) => c1 = next,
            None => return None,
        }
        //enter c1
        // this plane is crossed at tm.x, so this becomes the t for next sub_voxel
        // crossed this plane, so we add the dt.x (cuz the mid plane of that voxel is at that t)
        if let Some(hitfo) = self.try_hit_subvoxel(c1, ray, tm.x, planes.0, planes, &dt, &Vec3d::new(tm.x + dt.x, tm.y - dt.y, tm.z - dt.z)) {
            return Some(hitfo)
        }

        // crossing the second mid plane
        let c2: u8;
        match self.get_next_voxel(c1, planes.1, planes) {
            Some(next) => c2 = next,
            None => return None,
        }
        //enter c2
        if let Some(hitfo) = self.try_hit_subvoxel(c2, ray, tm.y, planes.1, planes, &dt, &Vec3d::new(tm.x + dt.x, tm.y + dt.y, tm.z - dt.z)) {
            return Some(hitfo)
        }

        // crossing the third plane
        let c3: u8;
        match self.get_next_voxel(c2, planes.2, planes) {
            Some(next) => c3 = next,
            None => return None,
        }
        //enter c3
        if let Some(hitfo) = self.try_hit_subvoxel(c3, ray, tm.z, planes.2, planes, &dt, &Vec3d::new(tm.x + dt.x, tm.y + dt.y, tm.z + dt.z)) {
            return Some(hitfo)
        }
        None
    }

    #[inline(always)]
    fn try_hit_subvoxel(&self, child_mask: u8, ray: &Ray, t: f64, t_plane: u8, planes: (u8, u8, u8), dt: &Vec3d, tm: &Vec3d) -> Option<RayHitfo> {
        if child_mask & self.branch_mask > 0 {
            let child = self.get_branch(OctreePos::new(child_mask));
            if let Some(hitfo) = child.hit(ray, t, t_plane, planes, dt, tm) {
                return Some(hitfo)
            }
        } else if child_mask & self.child_mask > 0 {
            // its a voxel. so return
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

    #[inline(always)]
    fn get_next_voxel(&self, child: u8, plane_mask: u8, planes: (u8, u8, u8)) -> Option<u8> {
        // gives the voxel_map when crossing a plane (after crossing the yz plane, we either go to the right voxel, or exit)
        if child & plane_mask > 0 {return None} // child is already in that side (eg -> move right (but im already in right, ig ill exit))
        // we need to know here what plane the plane_mask refers to
        let mut other_masks = {
            if plane_mask == planes.0 || plane_mask == !planes.0 {(planes.1, planes.2)}
            else if plane_mask == planes.1 || plane_mask == !planes.1 {(planes.0, planes.2)}
            else {(planes.0, planes.1)}
        };
        other_masks.0 = if other_masks.0 & child > 0 {other_masks.0} else {!other_masks.0};
        other_masks.1 = if other_masks.1 & child > 0 {other_masks.1} else {!other_masks.1};
        let next = plane_mask & other_masks.0 & other_masks.1;
        // eg: to more right from LDB, we first find the mask of every voxel with DB and "&" it with plane_mask(R) to the the voxel towards right
        Some(next)
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

