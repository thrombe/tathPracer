
use super::vec3d::Vec3d;
use super::math;

use super::material::{Material, Lit};
use super::ray::{Ray, RayHitfo};
use super::world::Rng;
use super::voxel_octree::{VoxelOctree, VoxelOctreeBranch, OctreePos, BbHit};
use super::object_octree::{ObjectOctree, ObjectOctreeBranch};
use super::triangle_octree::{Triangle, TriangleOctree, TriangleOctreeBranch};
use super::aabb::Aabb;

#[derive(Debug)]
pub enum Object {
    Sphere(Sphere),
    Plane(Plane),
    VoxelOctree(VoxelOctree),
    ObjectOctree(ObjectOctree),
    TriangleOctree(TriangleOctree),
}

impl Object {
    #[inline(always)]
    pub fn hit(&self, ray: &Ray, t_correction: f64, rng: &mut Rng) -> Option<RayHitfo> {
        match self {
            Object::Sphere(obj) => obj.hit(ray, t_correction),
            Object::Plane(obj) => obj.hit(ray, t_correction),
            Object::VoxelOctree(obj) => obj.hit(ray, t_correction, rng),
            Object::ObjectOctree(obj) => obj.hit(ray, t_correction, rng),
            Object::TriangleOctree(obj) => obj.hit(ray, t_correction),
        }
    }

    // octree size only cuz the plane has to go inside the main tree? 
    pub fn get_aabb(&self) -> Aabb {
        match self {
            Object::Sphere(obj) => obj.get_aabb(),
            Object::Plane(obj) => obj.get_aabb(),
            Object::VoxelOctree(obj) => obj.get_aabb(),
            Object::ObjectOctree(obj) => obj.get_aabb(),
            Object::TriangleOctree(obj) => obj.get_aabb(),
        }
    }
}

#[derive(Debug)]
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

    pub fn get_aabb(&self) -> Aabb {
        let radius = self.radius.abs();
        let half_diag = Vec3d::new(radius, radius, radius);
        Aabb::new(self.center-half_diag, self.center+half_diag)
    }
}

#[derive(Debug)]
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

    #[inline(always)]
    pub fn get_aabb(&self) -> Aabb {
        Aabb::new(Vec3d::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY), Vec3d::new(f64::INFINITY, f64::INFINITY, f64::INFINITY))
    }
}

impl VoxelOctree {
    #[inline]
    pub fn hit(&self, ray: &Ray, t_correction: f64, rng: &mut Rng) -> Option<RayHitfo> {
        // intersect the Nabb and return if not hit
        // todo

        let ray = Ray::new(self.world_to_tree_space(ray.pos), ray.dir);
        let mut planes = OctreePos::get_rub_masks();
        let bbox_min = Vec3d::new(-1.0, -1.0, -1.0);
        let (dt, t0) = {
            let mut ray = ray.clone();
            // the algo only works for +ve ray.dir, so if -ve, flip stuff
            // flip the ray wrt the mid planes (local)
            // flip the plane (call the left as right and right as left)
            if ray.dir.x < 0.0 {
                ray.dir.x *= -1.0;
                ray.pos.x *= -1.0;
                planes.0 = !planes.0;
            }
            if ray.dir.y < 0.0 {
                ray.dir.y *= -1.0;
                ray.pos.y *= -1.0;
                planes.1 = !planes.1;
            }
            if ray.dir.z < 0.0 {
                ray.dir.z *= -1.0;
                ray.pos.z *= -1.0;
                planes.2 = !planes.2;
            }
            let dt = Vec3d::new(1.0/ray.dir.x, 1.0/ray.dir.y, 1.0/ray.dir.z); // f64::Infinity comes up if any is 0.0
            let mut t0 = Vec3d::new((bbox_min.x-ray.pos.x)*dt.x, (bbox_min.y-ray.pos.y)*dt.y, (bbox_min.z-ray.pos.z)*dt.z);
            // inf is okay, nut nan is a problem. here, the nan is cuz 0*inf, so we can have 0 instead since this means ray.pos is on the line (prev line)
            // -ve inf is also a prob. -ve inf + inf -> nan. turning -ve inf to inf prevents the nan and also works out for possible good rays cuz max(t0) == inf is never a good ray
            t0.x = if t0.x.is_nan() {0.0} else if t0.x == f64::NEG_INFINITY {f64::INFINITY} else {t0.x};
            t0.y = if t0.y.is_nan() {0.0} else if t0.y == f64::NEG_INFINITY {f64::INFINITY} else {t0.y};
            t0.z = if t0.z.is_nan() {0.0} else if t0.z == f64::NEG_INFINITY {f64::INFINITY} else {t0.z};

            let t0 = (BbHit::new(t0.x, planes.0), BbHit::new(t0.y, planes.1), BbHit::new(t0.z, planes.2));
            (dt, t0)
        };
        { // return if it dosent hit the octree or octree behind the ray origin
            let t1 = (t0.0+dt.x*2.0, t0.1+dt.y*2.0, t0.2+dt.z*2.0);
            let min_t1 = math::min_vec(vec![t1.0.t, t1.1.t, t1.2.t]);
            let max_t0 = math::max_vec(vec![t0.0.t, t0.1.t, t0.2.t]);
            if max_t0 >= min_t1 {return None}
            if min_t1 < 0.0 {return None} // this means both max_t0 and min_t1 are -ve, so the octree must be brhind the ray
        }

        let t = { // for first hit, t = max(t0)
            let mut t = t0.0;
            if t.t < t0.1.t {t = t0.1}
            if t.t < t0.2.t {t = t0.2}
            t
        };

        // {
        //     if t.t == f64::INFINITY {
        //         dbg!(t1, t0, t, dt, max_t0, min_t1, max_t0 > min_t1);
        //         panic!();
        //     }
        // }
        
        // t value for first voxel is max(t for bbox.min)
        // how? -> well, the ray is in the voxel only if its inside all 3 plane boundaries. so max(t0) ensures this
        if let Some(mut hitfo) = self.main_branch.hit(&ray, t, t0, &dt, t_correction, self.lod_depth_limit.map(|x| x+1)) { // +1 to get this depth and insert depth in same level
            if hitfo.2 { // if ray.pos inside transparent voxels, then shoot another ray till it hits some other material
                if let Some(hitfo1) = self.main_branch.hit_transparent(&ray, t, 0.0, t0, &dt, t_correction, self.lod_depth_limit.map(|x| x+1), &self.materials[hitfo.1 as usize], hitfo.1, rng) {
                    hitfo.0 = hitfo1.0;
                    hitfo.1 = hitfo1.1;
                }
            }

            Some(RayHitfo { // get new material, transform the ray back and send it off for more bounces
               ray: Ray::new(self.tree_to_world_space(hitfo.0.ray.pos), hitfo.0.ray.dir),
               material: self.materials[hitfo.1 as usize].clone(),
               ..hitfo.0
           })
        } else {
            None
        }
    }

    pub fn get_aabb(&self) -> Aabb {
        let half_diag = Vec3d::new(self.half_size, self.half_size, self.half_size);
        Aabb::new(self.center-half_diag, self.center+half_diag)
    }
}

impl VoxelOctreeBranch {
    // think about all explanation in this func as if ray.dir.x and y and z > 0
    // the planes are masks of voxels to right/up/down/whatever. so if we "&" multiple of these, we can constrict where the required voxel is
    // eg: let 000 represent x, y, z, so 111 & 010 would mean gimme something (out of 111) that has y = 1 (this makes more sense with u8 but im lazy)

    // dt is current branch's (t1-t0)/2
    // in the return type, the u16 is material_index, and the bool inficates whether the ray.pos is inside a transparent voxel or not
    pub fn hit(&self, ray: &Ray, t: BbHit, t0: (BbHit, BbHit, BbHit), dt: &Vec3d, t_correction: f64, mut depth: Option<u16>) -> Option<(RayHitfo, u16, bool)> {
        // dbg!(self.pos);

        depth = depth.map(|x| x-1);        

        // since the planes are parallel and seperated by a constant amount, we just need some addition
        // to figure out the next t-values
        let tm = (t0.0 + dt.x, t0.1 + dt.y, t0.2 + dt.z); // -ve inf + inf = nan -> set to -ve inf or inf?
        
        // the next t-values can be found just* by sorting the ts in ascending order, since the ray hits the voxel in that order 
        let t1 = (tm.0 + dt.x, tm.1 + dt.y, tm.2 + dt.z);
        let ts = t.get_next_hits(vec![tm.0, tm.1, tm.2, t1.0, t1.1, t1.2]);
        
        let dt_by_2 = *dt*0.5; // dt for every sub voxel is half of its parents (since size is half)
        
        let entry_child_mask = if tm.0.t < t.t {tm.0.plane} else {!tm.0.plane}
                             & if tm.1.t < t.t {tm.1.plane} else {!tm.1.plane}
                             & if tm.2.t < t.t {tm.2.plane} else {!tm.2.plane};
        
        // now try and visit every child in this branch which could be hit (in order of ray intersection)
        // at most 4 sub-voxels can be hit
        
        if let Some(hitfo) = self.try_hit_subvoxel(entry_child_mask, ray, t, ts[0], t0, dt, &dt_by_2, t_correction, depth) {
            return Some(hitfo)
        }
        
        let mut child = entry_child_mask;
        for i in 0..3 {
            match self.get_next_voxel(child, ts[i], t0) {
                Some(next) => child = next,
                None => return None, // the ray exited this voxel(this branch)
            }
            if let Some(hitfo) = self.try_hit_subvoxel(child, ray, ts[i], ts[i+1], t0, dt, &dt_by_2, t_correction, depth) {
                return Some(hitfo)
            }
        }
        None
    }

    #[inline(always)]
    pub fn try_hit_subvoxel(&self, child_mask: u8, ray: &Ray, t: BbHit, ts_p1: BbHit, t0: (BbHit, BbHit, BbHit), dt: &Vec3d, dt_by_2: &Vec3d, t_correction: f64, depth: Option<u16>) -> Option<(RayHitfo, u16, bool)> {
        if ts_p1.t < t_correction {return None} // since ts > t and ts < 0 -> this voxel is completely behind the ray, so no need to enter
        if depth != Some(0) && child_mask & self.branch_mask > 0 { // check if the voxel is a branch
            let child = self.get_branch(child_mask);
            let t0 = self.get_t0_for(child_mask, dt, t0);
            if let Some(hitfo) = child.hit(ray, t, t0, dt_by_2, t_correction, depth) {
                return Some(hitfo)
            }
        } else if child_mask & self.child_mask > 0 { // check if leaf
            if t.t < t_correction { // if ray originates from somewhere in octree, we need to ignore -ve t. but we cant ignore the non-leaves if t -ve for them
                if self.transparency_mask & child_mask > 0 && ts_p1.t > t_correction { // this means ray originates from inside a transparent voxel thing
                    let hitfo = self.submit_hitfo(ray, BbHit::new(0.0, t.plane), child_mask, false, self.get_material_index(child_mask));
                    return Some((hitfo.0, hitfo.1, true))
                } else {
                    return None
                }    
            }
            let hitfo = self.submit_hitfo(ray, t, child_mask, false, self.get_material_index(child_mask));
            return Some((hitfo.0, hitfo.1, false))
        }
        None
    }

    // same as normal hit func but some differences for transparent voxels
    pub fn hit_transparent(&self, ray: &Ray, t: BbHit, mut min_t1: f64, t0: (BbHit, BbHit, BbHit), dt: &Vec3d, t_correction: f64, mut depth: Option<u16>, current_material: &Material, current_material_index: u16, rng: &mut Rng) -> Option<(RayHitfo, u16)> {
        depth = depth.map(|x| x-1);
        let tm = (t0.0 + dt.x, t0.1 + dt.y, t0.2 + dt.z);
        let t1 = (tm.0 + dt.x, tm.1 + dt.y, tm.2 + dt.z);
        match self.pos { // this is used to find the last voxel that the ray might encounter before it exits so that it can hit it
            OctreePos::Main => min_t1 = math::min_vec(vec![t1.0.t, t1.1.t, t1.2.t]),
            _ => (),
        }
        let ts = t.get_next_hits(vec![tm.0, tm.1, tm.2, t1.0, t1.1, t1.2]);
        let dt_by_2 = *dt*0.5;

        let entry_child_mask = if tm.0.t < t.t {tm.0.plane} else {!tm.0.plane}
                             & if tm.1.t < t.t {tm.1.plane} else {!tm.1.plane}
                             & if tm.2.t < t.t {tm.2.plane} else {!tm.2.plane};
        if let Some(hitfo) = self.try_hit_transparent_subvoxel(entry_child_mask, ray, t, ts[0], min_t1, t0, dt, &dt_by_2, t_correction, depth, current_material, current_material_index, rng) {
            return Some(hitfo)
        }
        
        let mut child = entry_child_mask;
        for i in 0..3 {
            match self.get_next_voxel(child, ts[i], t0) {
                Some(next) => child = next,
                None => return None,
            }
            if let Some(hitfo) = self.try_hit_transparent_subvoxel(child, ray, ts[i], ts[i+1], min_t1, t0, dt, &dt_by_2, t_correction, depth, current_material, current_material_index, rng) {
                return Some(hitfo)
            }
        }
        None
    }

    // this func ignores the voxel that are the same type as the voxel which has ray.pos
    #[inline(always)]
    pub fn try_hit_transparent_subvoxel(&self, child_mask: u8, ray: &Ray, t: BbHit, ts_p1: BbHit, min_t1: f64, t0: (BbHit, BbHit, BbHit), dt: &Vec3d, dt_by_2: &Vec3d, t_correction: f64, depth: Option<u16>, current_material: &Material, current_material_index: u16, rng: &mut Rng) -> Option<(RayHitfo, u16)> {
        if ts_p1.t < t_correction {return None}
        if depth != Some(0) && child_mask & self.branch_mask > 0 {
            let child = self.get_branch(child_mask);
            let t0 = self.get_t0_for(child_mask, dt, t0);
            if let Some(hitfo) = child.hit_transparent(ray, t, min_t1, t0, dt_by_2, t_correction, depth, current_material, current_material_index, rng) {
                return Some(hitfo)
            }
        } else if child_mask & self.child_mask > 0 { // if its a filled voxel, then it can be either the same type as the voxel with ray.pos or some other type
            if self.get_material_index(child_mask) == current_material_index {
                if (ts_p1.t - min_t1).abs() <= t_correction { // if the ray exits the entire octree, pretend an air voxel there
                    return Some(self.submit_hitfo(ray, ts_p1, child_mask, true, current_material_index))
                } else {
                    return None // if this voxel has the same material as the one with ray.pos, then ignore it
                }
            }

            // filled voxel which is not current material
            let mut t = t;
            let mut ray = ray.clone();
            ray.new_pos(t.t);
            let normal = {
                if self.normal_mask & child_mask > 0 {
                    self.get_normal(child_mask).clone()
                } else {
                    t.get_plane_normal()
                }
            };

            // if the next one is also transparent but of diff type, 2 refractions are needed
            if self.transparency_mask & child_mask > 0 {
                match current_material { // only transparent materials handled here
                    Material::Dielectric(_) => (),
                    _ => panic!("material not good"),
                }
                if let Some(rayy) = current_material.scatter(&ray, &(normal*(-1.0)), rng) {
                    ray = rayy;
                } else {panic!("ray not found")}
                t.t += t_correction; // to make sure its inside the new transparent voxel
            } else {
                t.t -= t_correction; // to make sure its not inside the non_transparent voxel. this prevents extra refractions
            }

            return Some((RayHitfo {
                t: t.t,
                normal, ray,
                material: Material::Lit(Lit {color: Vec3d::zero()}),
            }, self.get_material_index(child_mask)))
    
        } else { // air blocks inside octree -> return hit with flipped normal
            return Some(self.submit_hitfo(ray, t, child_mask, true, current_material_index))
        }
        None
    }

    #[inline(always)]
    pub fn submit_hitfo(&self, ray: &Ray, t: BbHit, child_mask: u8, flip_normal: bool, material_index: u16) -> (RayHitfo, u16) {
        let mut ray = ray.clone();
        ray.new_pos(t.t);
        let mut normal = {
            if self.normal_mask & child_mask > 0 {
                self.get_normal(child_mask).clone()
            } else {
                t.get_plane_normal()
            }
        };
        if flip_normal {normal *= -1.0}
        (RayHitfo {
            t: t.t,
            normal,
            material: Material::Lit(Lit {color: Vec3d::zero()}),
            ray,
        }, material_index)
    }

    /// gives the voxel_pos when crossing a plane (after crossing the yz plane, we either go to the right voxel, or exit)
    /// t is only used fot the planes. and not for the t values
    /// t_passed is the latest t value by which the ray passed a plane
    #[inline(always)]
    pub fn get_next_voxel(&self, child: u8, t_passed: BbHit, t: (BbHit, BbHit, BbHit)) -> Option<u8> {
        // eg: to more right from LDB, we first find the mask of every voxel with DB and "&" it with plane_mask(R) to the the voxel towards right
        let planes = (t.0.plane, t.1.plane, t.2.plane);
        let plane_mask = t_passed.plane;
        if child & plane_mask > 0 {return None} // child is already in that side (eg -> move right (but im already in right, ig ill exit))
        let mut other_masks = { // find the masks for planes other than t_passed
            if plane_mask == planes.0 || plane_mask == !planes.0 {(planes.1, planes.2)}
            else if plane_mask == planes.1 || plane_mask == !planes.1 {(planes.0, planes.2)}
            else {(planes.0, planes.1)}
        };
        other_masks.0 = if other_masks.0 & child > 0 {other_masks.0} else {!other_masks.0};
        other_masks.1 = if other_masks.1 & child > 0 {other_masks.1} else {!other_masks.1};
        let next = plane_mask & other_masks.0 & other_masks.1;
        Some(next)
    }

    // dt is of parent branch, not child
    #[inline(always)]
    pub fn get_t0_for(&self, child: u8, dt: &Vec3d, t0: (BbHit, BbHit, BbHit)) -> (BbHit, BbHit, BbHit) {
        // the t0 of branch can be reused here, since the child bbox_min coords only differ in 3 axes at most
        // and the bbox planes are only offset from the branch by 0 or side/2
        let mut t = t0;
        // if the voxel lies somewhere if the right side, we know that the t0.0 (yz plane) will be offset by side/2
        if child & t.0.plane > 0 { // right
            t.0.t += dt.x;
        } else { // left
        }
        if child & t.1.plane > 0 { // up
            t.1.t += dt.y;
        } else { // down
        }
        if child & t.2.plane > 0 { // back
            t.2.t += dt.z;
        } else { // fwd
        }
        t
    }
}

impl ObjectOctree {

    // look in voxel_octree for explanation of algorithm
    // only the differences are commented
    #[inline]
    pub fn hit(&self, ray_world_space: &Ray, t_correction: f64, rng: &mut Rng) -> Option<RayHitfo> {

        let ray = Ray::new(self.world_to_tree_space(ray_world_space.pos), ray_world_space.dir);
        let mut planes = OctreePos::get_rub_masks();
        let bbox_min = Vec3d::new(-1.0, -1.0, -1.0);
        let (dt, t0) = {
            let mut ray = ray.clone();
            if ray.dir.x < 0.0 {
                ray.dir.x *= -1.0;
                ray.pos.x *= -1.0;
                planes.0 = !planes.0;
            }
            if ray.dir.y < 0.0 {
                ray.dir.y *= -1.0;
                ray.pos.y *= -1.0;
                planes.1 = !planes.1;
            }
            if ray.dir.z < 0.0 {
                ray.dir.z *= -1.0;
                ray.pos.z *= -1.0;
                planes.2 = !planes.2;
            }
            let dt = Vec3d::new(1.0/ray.dir.x, 1.0/ray.dir.y, 1.0/ray.dir.z);
            let mut t0 = Vec3d::new((bbox_min.x-ray.pos.x)*dt.x, (bbox_min.y-ray.pos.y)*dt.y, (bbox_min.z-ray.pos.z)*dt.z);
            t0.x = if t0.x.is_nan() {0.0} else if t0.x == f64::NEG_INFINITY {f64::INFINITY} else {t0.x};
            t0.y = if t0.y.is_nan() {0.0} else if t0.y == f64::NEG_INFINITY {f64::INFINITY} else {t0.y};
            t0.z = if t0.z.is_nan() {0.0} else if t0.z == f64::NEG_INFINITY {f64::INFINITY} else {t0.z};

            let t0 = (BbHit::new(t0.x, planes.0), BbHit::new(t0.y, planes.1), BbHit::new(t0.z, planes.2));
            (dt, t0)
        };
        { // return if it dosent hit the octree or octree behind the ray origin
            let t1 = (t0.0+dt.x*2.0, t0.1+dt.y*2.0, t0.2+dt.z*2.0);
            let min_t1 = math::min_vec(vec![t1.0.t, t1.1.t, t1.2.t]);
            let max_t0 = math::max_vec(vec![t0.0.t, t0.1.t, t0.2.t]);
            if max_t0 >= min_t1 {return None}
            if min_t1 < 0.0 {return None}
        }

        let t = {
            let mut t = t0.0;
            if t.t < t0.1.t {t = t0.1}
            if t.t < t0.2.t {t = t0.2}
            t
        };

        // here sending ray_world_space cuz the hit funcs of objects expect them to be in world space, not tree space
        self.main_branch.hit(&ray_world_space, t, t0, &dt, t_correction, rng)
    }

    pub fn get_aabb(&self) -> Aabb {
        let half_diag = Vec3d::new(self.half_size, self.half_size, self.half_size);
        Aabb::new(self.center-half_diag, self.center+half_diag)
    }
}

impl ObjectOctreeBranch {
    
    // look in voxel_octree for explanation of algorithm
    pub fn hit(&self, ray: &Ray, t: BbHit, t0: (BbHit, BbHit, BbHit), dt: &Vec3d, t_correction: f64, rng: &mut Rng) -> Option<RayHitfo> {        
        let tm = (t0.0 + dt.x, t0.1 + dt.y, t0.2 + dt.z);
        let t1 = (tm.0 + dt.x, tm.1 + dt.y, tm.2 + dt.z);
        let ts = t.get_next_hits(vec![tm.0, tm.1, tm.2, t1.0, t1.1, t1.2]);
        let dt_by_2 = *dt*0.5;
        let entry_child_mask = if tm.0.t < t.t {tm.0.plane} else {!tm.0.plane}
                             & if tm.1.t < t.t {tm.1.plane} else {!tm.1.plane}
                             & if tm.2.t < t.t {tm.2.plane} else {!tm.2.plane};
        
        // to keep track of the closest hit within this voxel
        let mut closest = None;
        let mut min_t = f64::INFINITY;
        let mut no_further = false; // if one of the subvoxels have something hit, then the next subvoxels are just gonna have a greater t. so skip these
        
        if let Some(hitfo) = self.try_hit_subvoxel(entry_child_mask, ray, t, ts[0], t0, dt, &dt_by_2, t_correction, rng) {
            min_t = hitfo.t;
            closest = Some(hitfo);
            no_further = true;
        }
        
        let mut child = entry_child_mask;
        for i in 0..3 {
            if no_further {break}
            match self.get_next_voxel(child, ts[i], t0) {
                Some(next) => child = next,
                None => break, // got out of current branch
            }
            if let Some(hitfo) = self.try_hit_subvoxel(child, ray, ts[i], ts[i+1], t0, dt, &dt_by_2, t_correction, rng) {
                if hitfo.t < min_t {
                    min_t = hitfo.t;
                    closest = Some(hitfo);        
                }
                no_further = true;
            }
        }

        // required t = min(objects in curr branch, min(objects in subvoxels))
        if let Some(hitfo) = self.hit_objects_inside(ray, t_correction, rng) {
            if hitfo.t < min_t {
                // min_t = hitfo.t;
                closest = Some(hitfo);        
            }
        }
        closest
    }

    #[inline(always)]
    pub fn try_hit_subvoxel(&self, child_mask: u8, ray: &Ray, t: BbHit, ts_p1: BbHit, t0: (BbHit, BbHit, BbHit), dt: &Vec3d, dt_by_2: &Vec3d, t_correction: f64, rng: &mut Rng) -> Option<RayHitfo> {
        if ts_p1.t < t_correction {return None}
        if child_mask & self.branch_mask > 0 { // check if the voxel is a branch
            let child = self.get_branch(child_mask);
            let t0 = self.get_t0_for(child_mask, dt, t0);
            if let Some(hitfo) = child.hit(ray, t, t0, dt_by_2, t_correction, rng) {
                return Some(hitfo)
            }
        }
        None
    }

    // find min(t of objects in this branch)
    #[inline(always)]
    fn hit_objects_inside(&self, ray: &Ray, t_correction: f64, rng: &mut Rng) -> Option<RayHitfo> {
        let mut min_t = f64::INFINITY;
        let mut closest = None;
        for obj in &self.objects {
            if let Some(hitfo) = obj.hit(ray, t_correction, rng) {
                if hitfo.t < min_t {
                    min_t = hitfo.t;
                    closest = Some(hitfo);
                }
            }
        }
        closest
    }

    // same as VoxelOctree
    #[inline(always)]
    pub fn get_next_voxel(&self, child: u8, t_passed: BbHit, t: (BbHit, BbHit, BbHit)) -> Option<u8> {
        let planes = (t.0.plane, t.1.plane, t.2.plane);
        let plane_mask = t_passed.plane;
        if child & plane_mask > 0 {return None}
        let mut other_masks = {
            if plane_mask == planes.0 || plane_mask == !planes.0 {(planes.1, planes.2)}
            else if plane_mask == planes.1 || plane_mask == !planes.1 {(planes.0, planes.2)}
            else {(planes.0, planes.1)}
        };
        other_masks.0 = if other_masks.0 & child > 0 {other_masks.0} else {!other_masks.0};
        other_masks.1 = if other_masks.1 & child > 0 {other_masks.1} else {!other_masks.1};
        let next = plane_mask & other_masks.0 & other_masks.1;
        Some(next)
    }

    // same as VoxelOctree
    #[inline(always)]
    pub fn get_t0_for(&self, child: u8, dt: &Vec3d, t0: (BbHit, BbHit, BbHit)) -> (BbHit, BbHit, BbHit) {
        let mut t = t0;
        if child & t.0.plane > 0 { // right
            t.0.t += dt.x;
        } else { // left
        }
        if child & t.1.plane > 0 { // up
            t.1.t += dt.y;
        } else { // down
        }
        if child & t.2.plane > 0 { // back
            t.2.t += dt.z;
        } else { // fwd
        }
        t
    }
}

impl Triangle {
    pub fn hit(&self, vertices: &Vec<Vec3d>, materials: &Vec<Material>, ray: &Ray, t_correction: f64) -> Option<RayHitfo> {
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
        let vertices = self.get_vertices(vertices);
        let t = (*vertices.0-ray.pos).dot(self.normal)/ray.dir.dot(self.normal);
        if t < t_correction {return None}

        let A = *vertices.0;
        let B = *vertices.1;
        let C = *vertices.2;
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

        let mut ray = ray.clone();
        ray.new_pos(t);
        Some(RayHitfo {
            ray, t,
            material: self.get_material(materials),
            normal: self.normal,
        })
    }

    pub fn get_aabb(&self, vertices: (&Vec3d, &Vec3d, &Vec3d)) -> Aabb {
        let min = Vec3d::new(
            math::min_vec(vec![vertices.0.x, vertices.1.x, vertices.2.x]),
            math::min_vec(vec![vertices.0.y, vertices.1.y, vertices.2.y]),
            math::min_vec(vec![vertices.0.z, vertices.1.z, vertices.2.z]),
        );
        let max = Vec3d::new(
            math::max_vec(vec![vertices.0.x, vertices.1.x, vertices.2.x]),
            math::max_vec(vec![vertices.0.y, vertices.1.y, vertices.2.y]),
            math::max_vec(vec![vertices.0.z, vertices.1.z, vertices.2.z]),
        );
        Aabb::new(min, max)
    }
}

impl TriangleOctree {
    // look in voxel_octree for explanation of algorithm
    // only the differences are commented
    #[inline]
    pub fn hit(&self, ray_world_space: &Ray, t_correction: f64) -> Option<RayHitfo> {

        let ray = Ray::new(self.world_to_tree_space(ray_world_space.pos), ray_world_space.dir);
        let mut planes = OctreePos::get_rub_masks();
        let bbox_min = Vec3d::new(-1.0, -1.0, -1.0);
        let (dt, t0) = {
            let mut ray = ray.clone();
            if ray.dir.x < 0.0 {
                ray.dir.x *= -1.0;
                ray.pos.x *= -1.0;
                planes.0 = !planes.0;
            }
            if ray.dir.y < 0.0 {
                ray.dir.y *= -1.0;
                ray.pos.y *= -1.0;
                planes.1 = !planes.1;
            }
            if ray.dir.z < 0.0 {
                ray.dir.z *= -1.0;
                ray.pos.z *= -1.0;
                planes.2 = !planes.2;
            }
            let dt = Vec3d::new(1.0/ray.dir.x, 1.0/ray.dir.y, 1.0/ray.dir.z);
            let mut t0 = Vec3d::new((bbox_min.x-ray.pos.x)*dt.x, (bbox_min.y-ray.pos.y)*dt.y, (bbox_min.z-ray.pos.z)*dt.z);
            t0.x = if t0.x.is_nan() {0.0} else if t0.x == f64::NEG_INFINITY {f64::INFINITY} else {t0.x};
            t0.y = if t0.y.is_nan() {0.0} else if t0.y == f64::NEG_INFINITY {f64::INFINITY} else {t0.y};
            t0.z = if t0.z.is_nan() {0.0} else if t0.z == f64::NEG_INFINITY {f64::INFINITY} else {t0.z};

            let t0 = (BbHit::new(t0.x, planes.0), BbHit::new(t0.y, planes.1), BbHit::new(t0.z, planes.2));
            (dt, t0)
        };
        { // return if it dosent hit the octree or octree behind the ray origin
            let t1 = (t0.0+dt.x*2.0, t0.1+dt.y*2.0, t0.2+dt.z*2.0);
            let min_t1 = math::min_vec(vec![t1.0.t, t1.1.t, t1.2.t]);
            let max_t0 = math::max_vec(vec![t0.0.t, t0.1.t, t0.2.t]);
            if max_t0 >= min_t1 {return None}
            if min_t1 < 0.0 {return None}
        }

        let t = {
            let mut t = t0.0;
            if t.t < t0.1.t {t = t0.1}
            if t.t < t0.2.t {t = t0.2}
            t
        };

        self.main_branch.hit(&ray_world_space, &self.vertices, &self.materials, t, t0, &dt, t_correction)
    }

    pub fn get_aabb(&self) -> Aabb {
        let half_diag = Vec3d::new(self.half_size, self.half_size, self.half_size);
        Aabb::new(self.center-half_diag, self.center+half_diag)
    }
}

impl TriangleOctreeBranch {
    // look in voxel_octree for explanation of algorithm
    pub fn hit(&self, ray: &Ray, vertices: &Vec<Vec3d>, materials: &Vec<Material>, t: BbHit, t0: (BbHit, BbHit, BbHit), dt: &Vec3d, t_correction: f64) -> Option<RayHitfo> {        
        let tm = (t0.0 + dt.x, t0.1 + dt.y, t0.2 + dt.z);
        let t1 = (tm.0 + dt.x, tm.1 + dt.y, tm.2 + dt.z);
        let ts = t.get_next_hits(vec![tm.0, tm.1, tm.2, t1.0, t1.1, t1.2]);
        let dt_by_2 = *dt*0.5;
        let entry_child_mask = if tm.0.t < t.t {tm.0.plane} else {!tm.0.plane}
                             & if tm.1.t < t.t {tm.1.plane} else {!tm.1.plane}
                             & if tm.2.t < t.t {tm.2.plane} else {!tm.2.plane};
        
        // to keep track of the closest hit within this voxel
        let mut closest = None;
        let mut min_t = f64::INFINITY;
        let mut no_further = false; // if one of the subvoxels have something hit, then the next subvoxels are just gonna have a greater t. so skip these
        
        if let Some(hitfo) = self.try_hit_subvoxel(entry_child_mask, ray, vertices, materials, t, ts[0], t0, dt, &dt_by_2, t_correction) {
            min_t = hitfo.t;
            closest = Some(hitfo);
            no_further = true;
        }
        
        let mut child = entry_child_mask;
        for i in 0..3 {
            if no_further {break}
            match self.get_next_voxel(child, ts[i], t0) {
                Some(next) => child = next,
                None => break, // got out of current branch
            }
            if let Some(hitfo) = self.try_hit_subvoxel(child, ray, vertices, materials, ts[i], ts[i+1], t0, dt, &dt_by_2, t_correction) {
                if hitfo.t < min_t {
                    min_t = hitfo.t;
                    closest = Some(hitfo);        
                }
                no_further = true;
            }
        }

        // required t = min(objects in curr branch, min(objects in subvoxels))
        if let Some(hitfo) = self.hit_objects_inside(ray, vertices, materials, t_correction) {
            if hitfo.t < min_t {
                // min_t = hitfo.t;
                closest = Some(hitfo);        
            }
        }
        closest
    }

    #[inline(always)]
    pub fn try_hit_subvoxel(&self, child_mask: u8, ray: &Ray, vertices: &Vec<Vec3d>, materials: &Vec<Material>, t: BbHit, ts_p1: BbHit, t0: (BbHit, BbHit, BbHit), dt: &Vec3d, dt_by_2: &Vec3d, t_correction: f64) -> Option<RayHitfo> {
        if ts_p1.t < t_correction {return None}
        if child_mask & self.branch_mask > 0 { // check if the voxel is a branch
            let child = self.get_branch(child_mask);
            let t0 = self.get_t0_for(child_mask, dt, t0);
            if let Some(hitfo) = child.hit(ray, vertices, materials, t, t0, dt_by_2, t_correction) {
                return Some(hitfo)
            }
        }
        None
    }

    // find min(t of objects in this branch)
    #[inline(always)]
    fn hit_objects_inside(&self, ray: &Ray, vertices: &Vec<Vec3d>, materials: &Vec<Material>, t_correction: f64) -> Option<RayHitfo> {
        let mut min_t = f64::INFINITY;
        let mut closest = None;
        for tri in &self.triangles {
            if let Some(hitfo) = tri.hit(vertices, materials, ray, t_correction) {
                if hitfo.t < min_t {
                    min_t = hitfo.t;
                    closest = Some(hitfo);
                }
            }
        }
        closest
    }

    // same as VoxelOctree
    #[inline(always)]
    pub fn get_next_voxel(&self, child: u8, t_passed: BbHit, t: (BbHit, BbHit, BbHit)) -> Option<u8> {
        let planes = (t.0.plane, t.1.plane, t.2.plane);
        let plane_mask = t_passed.plane;
        if child & plane_mask > 0 {return None}
        let mut other_masks = {
            if plane_mask == planes.0 || plane_mask == !planes.0 {(planes.1, planes.2)}
            else if plane_mask == planes.1 || plane_mask == !planes.1 {(planes.0, planes.2)}
            else {(planes.0, planes.1)}
        };
        other_masks.0 = if other_masks.0 & child > 0 {other_masks.0} else {!other_masks.0};
        other_masks.1 = if other_masks.1 & child > 0 {other_masks.1} else {!other_masks.1};
        let next = plane_mask & other_masks.0 & other_masks.1;
        Some(next)
    }

    // same as VoxelOctree
    #[inline(always)]
    pub fn get_t0_for(&self, child: u8, dt: &Vec3d, t0: (BbHit, BbHit, BbHit)) -> (BbHit, BbHit, BbHit) {
        let mut t = t0;
        if child & t.0.plane > 0 { // right
            t.0.t += dt.x;
        } else { // left
        }
        if child & t.1.plane > 0 { // up
            t.1.t += dt.y;
        } else { // down
        }
        if child & t.2.plane > 0 { // back
            t.2.t += dt.z;
        } else { // fwd
        }
        t
    }
}