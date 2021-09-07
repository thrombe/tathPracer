
use super::vec3d::Vec3d;
use super::math;
use super::voxel_octree::{OctreePos, BbHit};
use super::objects::Object;
use super::aabb::Aabb;
use super::ray::{Ray, RayHitfo};
use super::world::Rng;

#[derive(Debug)]
pub struct WorldOctree {
    // entire octree lies in -1 to 1, half_size converts to and from this to world space or whatever
    half_size: f64,
    pub scale_factor: f64,
    pub main_branch: WorldOctreeBranch,
    // ! object octree plz
}

impl WorldOctree {
    pub fn new(size: f64) -> Self {
        Self {
            half_size: size/2.0,
            scale_factor: 2.0/size,
            main_branch: WorldOctreeBranch::new(OctreePos::Main),
        }
    }

    pub fn insert_object(&mut self, object: Object) {
        let aabb = object.get_aabb();
        if aabb.min.x == f64::NEG_INFINITY && aabb.max.y == f64::INFINITY { // for planes and stuff
            self.main_branch.objects.push(object);
            return
        }
        let tree_aabb = Aabb::new(Vec3d::new(-self.half_size, -self.half_size, -self.half_size), Vec3d::new(self.half_size, self.half_size, self.half_size));
        if !(tree_aabb.contains(&aabb)) {
            dbg!(&object);
            panic!("object out of tree");
        }
        self.main_branch.insert_object(object, aabb*self.scale_factor);
    }

    #[inline(always)]
    pub fn world_to_tree_space(&self, pos: Vec3d) -> Vec3d {
        pos*self.scale_factor
    }

    #[inline(always)]
    pub fn tree_to_world_space(&self, pos: Vec3d) -> Vec3d {
        pos*self.half_size
    }

    /// deletes empty branches, and stuff
    pub fn compress_sparseness(&mut self) {
        todo!();
    }
}

#[derive(Debug)]
pub struct WorldOctreeBranch {
    // indices of everything can be accessed the same way as VoxelOctree
    pub pos: OctreePos,
    pub child_mask: u8,
    pub branch_mask: u8,
    pub chilranches: Vec<WorldOctreeBranch>,
    pub objects: Vec<Object>,
}

impl WorldOctreeBranch { // all branches consider their space as -1 to 1
    fn new(pos: OctreePos) -> Self {
        Self {
            pos, child_mask: 0, branch_mask: 0, chilranches: Vec::<Self>::new(), objects: Vec::<Object>::new(),
        }
    }

    // bias towards +ve direction(if on line)
    fn get_pos_from_point(&self, point: &Vec3d) -> OctreePos {
        // println!("point {:?}", point);
        let (r, u, b) = OctreePos::get_rub_masks();
        let mut pos = if point.x >= 0.0 {r} else {!r};
        pos &= if point.y >= 0.0 {u} else {!u};
        pos &= if point.z >= 0.0 {b} else {!b};
        OctreePos::new(pos)
    }

    // center of subbranch in self space
    #[inline(always)]
    fn get_branch_coord(&self, pos: OctreePos) -> Vec3d {
        match pos { // r, u, b -> +ve, l, d, f -> -ve
            OctreePos::RUB => Vec3d::new(0.5, 0.5, 0.5),
            OctreePos::LUB => Vec3d::new(-0.5, 0.5, 0.5),
            OctreePos::LUF => Vec3d::new(-0.5, 0.5, -0.5),
            OctreePos::RUF => Vec3d::new(0.5, 0.5, -0.5),
            OctreePos::RDF => Vec3d::new(0.5, -0.5, -0.5),
            OctreePos::LDF => Vec3d::new(-0.5, -0.5, -0.5),
            OctreePos::LDB => Vec3d::new(-0.5, -0.5, 0.5),
            OctreePos::RDB => Vec3d::new(0.5, -0.5, 0.5),
            OctreePos::Main => Vec3d::zero(),
        }
    }

    // in parent's space
    fn get_coord(&self) -> Vec3d {
        self.get_branch_coord(self.pos)
    }

    #[inline(always)]
    fn self_aabb(&self) -> Aabb {
        Aabb::new(Vec3d::new(-1.0, -1.0, -1.0), Vec3d::new(1.0, 1.0, 1.0))
    }
    
    #[inline(always)]
    fn child_aabb(&self, pos: OctreePos) -> Aabb {
        let p5 = Vec3d::new(0.5, 0.5, 0.5);
        let c = self.get_branch_coord(pos);
        Aabb::new(c-p5, c+p5)
    }

    fn insert_object(&mut self, object: Object, aabb: Aabb) {
        let mut child_mask = 1;
        for _ in 0..8 {
            let pos = OctreePos::new(child_mask);
            if self.child_aabb(pos).contains(&aabb) {
                let aabb2 = (aabb-self.get_branch_coord(pos))*2.0;
                let branch = self.add_or_get_branch(child_mask);
                branch.insert_object(object, aabb2);
                return
            }
            child_mask *= 2;
        }
        self.objects.push(object);
    }

    fn add_or_get_branch(&mut self, pos_mask: u8) -> &mut Self {
        if self.child_mask & pos_mask > 0 {
            self.get_branch_mut(pos_mask)
        } else {
            self.add_branch(pos_mask)
        }
    }

    fn add_branch(&mut self, pos_mask: u8) -> &mut Self {
        // create a new voxel
        self.child_mask = self.child_mask | pos_mask;

        // create a new branch
        self.branch_mask = self.branch_mask | pos_mask;
        let chilranch = Self::new(OctreePos::new(pos_mask));

        // insert branch at correct index
        let index = self.get_info_index(self.branch_mask, pos_mask);
        self.chilranches.insert(index, chilranch);
        &mut self.chilranches[index]
    }

    #[inline(always)]
    pub fn get_branch_mut(&mut self, pos_mask: u8) -> &mut Self {
        let index = self.get_info_index(self.branch_mask, pos_mask);
        &mut self.chilranches[index]
    }
    
    #[inline(always)]
    pub fn get_branch(&self, pos_mask: u8) -> &Self {
        let index = self.get_info_index(self.branch_mask, pos_mask);
        &self.chilranches[index]
    }

    /// assuming that there is info for this voxel
    #[inline(always)]
    pub fn get_info_index(&self, info_mask: u8, pos_mask: u8) -> usize {
        let mut count = 0;
        let mut mask = 1;
        for _ in 0..8 {
            if mask == pos_mask {return count}
            if (info_mask & mask) > 0 {count += 1}
            mask *= 2; // move the on bit left
        }
        count // code never gets here, but rust complains. smh my head
    }
}

// hit functions below

impl WorldOctree {

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
}

impl WorldOctreeBranch {
    
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

    // dt is of parent branch, not child
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