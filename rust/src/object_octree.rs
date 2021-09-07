
use super::vec3d::Vec3d;
use super::voxel_octree::{OctreePos};
use super::objects::Object;
use super::aabb::Aabb;

#[derive(Debug)]
pub struct ObjectOctree {
    // entire octree lies in -1 to 1, half_size converts to and from this to world space or whatever
    pub half_size: f64,
    pub scale_factor: f64,
    pub center: Vec3d,
    pub main_branch: ObjectOctreeBranch,
}

impl ObjectOctree {
    pub fn new(center:Vec3d, size: f64) -> Self {
        Self {
            half_size: size/2.0,
            scale_factor: 2.0/size,
            center,
            main_branch: ObjectOctreeBranch::new(OctreePos::Main),
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
        (pos-self.center)*self.scale_factor
    }

    #[inline(always)]
    pub fn tree_to_world_space(&self, pos: Vec3d) -> Vec3d {
        (pos+self.center)*self.half_size
    }

    /// deletes empty branches, and stuff
    pub fn compress_sparseness(&mut self) {
        todo!();
    }
}

#[derive(Debug)]
pub struct ObjectOctreeBranch {
    // indices of everything can be accessed the same way as VoxelOctree
    pub pos: OctreePos,
    pub child_mask: u8,
    pub branch_mask: u8,
    pub chilranches: Vec<ObjectOctreeBranch>,
    pub objects: Vec<Object>,
}

impl ObjectOctreeBranch { // all branches consider their space as -1 to 1
    fn new(pos: OctreePos) -> Self {
        Self {
            pos, child_mask: 0, branch_mask: 0, chilranches: Vec::<Self>::new(), objects: Vec::<Object>::new(),
        }
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
}

// the methods which are exact copy of voxel_octree
impl ObjectOctreeBranch {

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