
// use std::sync::{Arc, RwLock};

use super::vec3d::Vec3d;
use super::material::{Material, Lit};
use super::math;

#[derive(Debug)]
pub struct Octree {
    // entire octree lies in -1 to 1, half_size converts to and from this to world space or whatever
    half_size: f64,
    pub scale_factor: f64,
    pub main_branch: OctreeBranch,
    pub materials: Vec<Material>, // material at 0 index should be some default material for undefined stuff
}

impl Octree {
    pub fn new(size: f64) -> Self {
        Self {
            half_size: size/2.0,
            scale_factor: 2.0/size,
            main_branch: OctreeBranch::new(OctreePos::Main),
            materials: vec![Material::Lit(Lit {color: Vec3d::zero()})],
        }
    }

    /// add new material using this, and use the returned index to fill in voxels
    pub fn add_material(&mut self, material: Material) -> u32 {
        self.materials.push(material);
        (self.materials.len()-1) as u32
    }

    pub fn insert_voxel(&mut self, mut point: Vec3d, depth: usize, material_index: u32, normal: Option<Vec3d>) {
        point = self.world_to_tree_space(point);
        if (math::abs(point.x) > 1.0) || (math::abs(point.y) > 1.0) || (math::abs(point.z) > 1.0) {panic!()}
        self.main_branch.insert_voxel_from_point(point, depth, material_index, normal);
    }

    pub fn world_to_tree_space(&self, pos: Vec3d) -> Vec3d {
        pos*self.scale_factor
    }

    pub fn tree_to_world_space(&self, pos: Vec3d) -> Vec3d {
        pos*self.half_size
    }

    /// deletes empty branches, deletes unnecessary voxels
    pub fn compress_sparseness(&mut self) {
        // check normals too!!
        todo!();
    }
}

#[derive(Debug)]
pub struct OctreeBranch {
    // children & branch_mask -> branched children, children & !branch_mask -> leafs, children -> live children
    pub pos: OctreePos,
    pub child_mask: u8,
    pub branch_mask: u8, // if branch, the bit is set
    pub chilranches: Vec<OctreeBranch>, // index -> same as normals but with branch_mask
    pub normal_mask: u8, // this space was being wasted anyway // if set, then this voxel has normal
    pub normals: Vec<Vec3d>, // index -> the index of the voxel in the normal_mask while ignoring 0's,
    // for eg, 01001000(normal_map), normal for 01000000(voxel) is 1 and for 00001000(voxel) is 0
    pub materials: Vec<u32>, // index from Octree.materials
}

impl OctreeBranch { // all branches consider their space as -1 to 1
    fn new(pos: OctreePos) -> Self {
        Self {
            pos, child_mask: 0, branch_mask: 0, chilranches: Vec::<Self>::new(), normal_mask: 0,
            normals: Vec::<Vec3d>::new(), materials: Vec::<u32>::new(),
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
    
    fn set_voxel(&mut self, pos: OctreePos, material_index: u32, normal: Option<Vec3d>) {
        let pos_mask = pos.get_mask();
        self.child_mask = self.child_mask | pos_mask;
        
        match normal {
            Some(normal) => self.insert_normal(normal, pos_mask),
            None => (),
        }
        
        let index = self.get_info_index(self.child_mask, pos_mask);
        self.materials.insert(index, material_index)
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

    fn insert_voxel_from_point(&mut self, point: Vec3d, depth: usize, material_index: u32, normal: Option<Vec3d>) {
        let pos = self.get_pos_from_point(&point);
        // println!("{:?}", pos);
        if depth == 0 {
            self.set_voxel(pos, material_index, normal);
            return
        }
        let branch = self.add_branch(pos, material_index, normal);
        branch.insert_voxel_from_point((point-branch.get_coord())*2.0, depth-1, material_index, normal);
    }
    
    fn add_branch(&mut self, pos: OctreePos, material_index: u32, normal: Option<Vec3d>) -> &mut Self {
        self.set_voxel(pos, material_index, normal);
        let pos_mask = pos.get_mask();
        if (self.branch_mask & pos_mask) > 0 {return self.get_branch_mut(pos)}
        self.branch_mask = self.branch_mask | pos_mask;
        let chilranch = Self::new(pos);
        // insert branch at correct index
        let index = self.get_info_index(self.branch_mask, pos_mask);
        self.chilranches.insert(index, chilranch);
        &mut self.chilranches[index]
    }

    // find where the normal should go (index), and correctly set/modify normals and indices
    fn insert_normal(&mut self, normal: Vec3d, pos_mask: u8) {
        self.normal_mask = self.normal_mask | pos_mask;
        let index = self.get_info_index(self.normal_mask, pos_mask);
        self.normals.insert(index, normal);
    }

    pub fn get_branch_mut(&mut self, pos: OctreePos) -> &mut Self {
        let index = self.get_info_index(self.branch_mask, pos.get_mask());
        &mut self.chilranches[index]
    }
    
    pub fn get_branch(&self, pos: OctreePos) -> &Self {
        let index = self.get_info_index(self.branch_mask, pos.get_mask());
        &self.chilranches[index]
    }

    pub fn get_normal(&self, pos_mask: u8) -> &Vec3d {
        &self.normals[self.get_info_index(self.normal_mask, pos_mask)]
    }

    /// assuming that there is info for this voxel
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

#[derive(Clone, Copy, Debug)]
pub struct BbHit { // bounding box hit
  pub t: f64,
  pub plane: u8,
}

impl BbHit {
    pub fn new(t: f64, plane: u8)-> Self {
        Self {t, plane}
    }

    // sort all ts in ascending order + remove any ts < current t
    pub fn get_next_hits(&self, mut hits: Vec<Self>) -> Vec<Self> {
        hits.sort_by(|a, b| a.t.partial_cmp(&b.t).unwrap());
        hits.iter().filter(|a| a.t > self.t).cloned().collect()
    }

    // same as above but only for the next emelent
    pub fn get_next_hit(&self, hits: Vec<Self>) -> Self {
        let mut best = BbHit::new(f64::INFINITY, 0);
        for hit in hits {
            if hit.t < best.t && hit.t > self.t {
                best = hit;
            }
        }
        if best.t == f64::INFINITY {panic!()}
        best
    }
}

use std::ops::{Add, Sub};
impl Add<f64> for BbHit {
    type Output = Self;

    #[inline(always)]
    fn add(self, t: f64) -> Self {
        Self::new(self.t + t, self.plane)
    }
}

impl Sub<f64> for BbHit {
    type Output = Self;

    #[inline(always)]
    fn sub(self, t: f64) -> Self {
        Self::new(self.t - t, self.plane)
    }
}
