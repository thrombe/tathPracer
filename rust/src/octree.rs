
use std::sync::{Arc, RwLock};

use super::vec3d::Vec3d;
use super::material::Material;
use super::math;

#[derive(Debug)]
pub struct Octree {
    // entire octree lies in -1 to 1, half_size converts to and from this to world space or whatever
    half_size: f64,
    pub scale_factor: f64,
    pub main_branch: OctreeBranch,
    // materials: Arc<Vec<Material>>, // do vecs need box to be in a struct?
    // materials: Vec<Material>,
}

impl Octree {
    pub fn new(size: f64) -> Self {
        Self {
            half_size: size/2.0,
            scale_factor: 2.0/size,
            main_branch: OctreeBranch::new(OctreePos::Main),
            // materials: Arc::new(Vec::<Material>::new()),
            // materials: Vec::<Material>::new(),
        }
    }

    pub fn insert_voxel(&mut self, mut point: Vec3d, depth: usize) {
        point = self.world_to_tree_space(point);
        if (math::abs(point.x) > 1.0) || (math::abs(point.y) > 1.0) || (math::abs(point.z) > 1.0) {panic!()}
        self.main_branch.insert_voxel_from_point(point, depth);
    }

    pub fn world_to_tree_space(&self, pos: Vec3d) -> Vec3d {
        pos*self.scale_factor
    }

    pub fn tree_to_world_space(&self, pos: Vec3d) -> Vec3d {
        pos*self.half_size
    }

    // pub fn tree_to_world_space_f64(&self, t: f64) -> f64 {
    //     t*self.half_size
    // } why is this needed?

    /// deletes empty branches, deletes unnecessary voxels
    pub fn compress_sparseness(&mut self) {
        // check normals too!!
    }
}

#[derive(Debug)]
pub struct OctreeBranch {
    // children & branch_mask -> branched children, children & !branch_mask -> leafs, children -> live children
    pub pos: OctreePos,
    pub child_mask: u8,
    pub branch_mask: u8, // if branch, the bit is set
    pub chilranches: Box<Vec<OctreeBranch>>, // index -> same as normals but with branch_mask
    // normal_mask: u8, // this space was being wasted anyway // if set, then this voxel has normal
    // normals: Arc<Vec<Vec3d>>, // index -> the index of the voxel in the normal_mask while ignoring 0's,
    // for eg, 01001000(normal_map), normal for 01000000(voxel) is 1 and for 00001000(voxel) is 0
    // material: u32, // index from Octree.materials
}

impl OctreeBranch { // all branches consider their space as -1 to 1
    fn new(pos: OctreePos) -> Self {
        Self {
            pos, child_mask: 0, branch_mask: 0, chilranches: Box::new(Vec::<Self>::new()),
        }
    }

    fn get_pos_from_point(&self, point: &Vec3d) -> OctreePos {
        // println!("point {:?}", point);
        let (r, u, b) = OctreePos::get_rub_masks();
        let mut pos = if point.x >= 0.0 {r} else {!r};
        pos &= if point.y >= 0.0 {u} else {!u};
        pos &= if point.z >= 0.0 {b} else {!b};
        OctreePos::new(pos)
    }
    
    fn set_voxel(&mut self, pos: OctreePos) {
        self.child_mask = self.child_mask | pos.get_mask();
    }

    // in self space
    fn get_voxel_coord(&self, pos: OctreePos) -> Vec3d {
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
        self.get_voxel_coord(self.pos)
    }

    fn insert_voxel_from_point(&mut self, point: Vec3d, depth: usize) {
        let pos = self.get_pos_from_point(&point);
        // println!("{:?}", pos);
        if depth == 0 {
            self.set_voxel(pos);
            return
        }
        let branch = self.add_branch(pos);
        branch.insert_voxel_from_point((point-branch.get_coord())*2.0, depth-1);
    }
    
    fn add_branch(&mut self, pos: OctreePos) -> &mut Self {
        self.set_voxel(pos);
        let pos_mask = pos.get_mask();
        if (self.branch_mask & pos_mask) > 0 {return self.get_branch_mut(pos)}
        self.branch_mask = self.branch_mask | pos_mask;
        let chilranch = Self::new(pos);
        // insert branch at correct index
        let index = self.get_branch_index(pos);
        self.chilranches.insert(index, chilranch);
        &mut self.chilranches[index]
    }

    pub fn get_branch_mut(&mut self, pos: OctreePos) -> &mut Self {
        let index = self.get_branch_index(pos);
        &mut self.chilranches[index]
    }

    pub fn get_branch(&self, pos: OctreePos) -> &Self {
        let index = self.get_branch_index(pos);
        &self.chilranches[index]
    }

    /// assuming that there is a branch for this
    fn get_branch_index(&self, pos: OctreePos) -> usize {
        let pos_mask = pos.get_mask();
        let mut count = 0;
        let mut mask = 1;
        for _ in 0..8 {
            if mask == pos_mask {return count}
            if (self.branch_mask & mask) > 0 {count += 1}
            mask *= 2; // move the on bit left
        }
        count // code never gets here, but rust complains. smh my head
    }
}

#[derive(Debug, Clone, Copy)]
pub enum OctreePos { // facing -z
    // U -> up, D -> down, F -> forward, B -> backward, L -> left, R -> right
    // 0 to 7 from left to right
    RUB, LUB, LUF, RUF, RDF, LDF, LDB, RDB, Main,
}

impl OctreePos {
    pub fn new(mask: u8) -> Self {
        match mask {
            0b_00000001 => OctreePos::RUB,
            0b_00000010 => OctreePos::LUB,
            0b_00000100 => OctreePos::LUF,
            0b_00001000 => OctreePos::RUF,
            0b_00010000 => OctreePos::RDF,
            0b_00100000 => OctreePos::LDF,
            0b_01000000 => OctreePos::LDB,
            0b_10000000 => OctreePos::RDB,
            _ => {
                println!("maskkkk --- {:?}", mask);
                panic!();
            },
        }
    }

    pub fn get_mask(&self) -> u8 {
        match self {
            OctreePos::RUB => 0b_00000001,
            OctreePos::LUB => 0b_00000010,
            OctreePos::LUF => 0b_00000100,
            OctreePos::RUF => 0b_00001000,
            OctreePos::RDF => 0b_00010000,
            OctreePos::LDF => 0b_00100000,
            OctreePos::LDB => 0b_01000000,
            OctreePos::RDB => 0b_10000000,
            OctreePos::Main => panic!(),
            // OctreePos::Exit => panic!(),
        }
    }
    
    pub fn get_rub_masks() -> (u8, u8, u8) {
        (
        0b_10011001, // R - 153 or ! 102
        0b_00001111, // U - 15 or ! 240
        0b_11000011, // B - 195 or ! 60
        )
    }
}

// #[derive(Debug)]
// struct Leaf<'a> {
//     pos: OctreePos,
//     normal: &'a Vec3d,
//     color: &'a Vec3d,
// }
