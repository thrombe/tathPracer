
// use std::sync::{Arc, RwLock};
// if multithreading gives trouble with arc and stuff, try the share/sync trait

use super::vec3d::Vec3d;
use super::material::{Material, Lit};
// use super::math;

#[derive(Debug)]
pub struct VoxelOctree {
    // entire octree lies in -1 to 1, half_size converts to and from this to world space or whatever
    pub half_size: f64,
    pub scale_factor: f64,
    pub center: Vec3d,
    pub main_branch: VoxelOctreeBranch,
    pub materials: Vec<Material>, // material at 0 index should be some default material for undefined stuff
    pub lod_depth_limit: Option<u16>,
}

impl VoxelOctree {
    pub fn new(center:Vec3d, size: f64, lod_depth_limit: Option<u16>) -> Self {
        Self {
            lod_depth_limit,
            center,
            half_size: size/2.0,
            scale_factor: 2.0/size,
            main_branch: VoxelOctreeBranch::new(OctreePos::Main),
            materials: vec![Material::Lit(Lit {color: Vec3d::zero()})],
        }
    }

    /// add new material using this, and use the returned index to fill in voxels
    pub fn add_material(&mut self, material: Material) -> u16 {
        self.materials.push(material);
        (self.materials.len()-1) as u16
    }

    pub fn insert_voxel(&mut self, mut point: Vec3d, depth: usize, material_index: u16, normal: Option<Vec3d>) {
        point = self.world_to_tree_space(point);
        let transparency = if let Material::Dielectric(_) = self.materials[material_index as usize] {true} else {false};
        if (point.x.abs() > 1.0) || (point.y.abs() > 1.0) || (point.z.abs() > 1.0) {panic!()}
        self.main_branch.insert_voxel_from_point(point, depth, material_index, normal, transparency);
    }

    #[inline(always)]
    pub fn world_to_tree_space(&self, pos: Vec3d) -> Vec3d {
        (pos-self.center)*self.scale_factor
    }

    #[inline(always)]
    pub fn tree_to_world_space(&self, pos: Vec3d) -> Vec3d {
        (pos+self.center)*self.half_size
    }

    /// deletes empty branches, deletes unnecessary voxels
    pub fn compress_sparseness(&mut self) {
        // check normals too!!
        // voxels with no normals and default materials can be considered mergable with voxels that have materials and normals (if not visible)
        todo!();
    }
}

#[derive(Debug)]
pub struct VoxelOctreeBranch {
    // children & branch_mask -> branched children, children & !branch_mask -> leafs, children -> live children
    pub pos: OctreePos,
    pub child_mask: u8,
    pub branch_mask: u8, // if branch, the bit is set
    pub chilranches: Vec<VoxelOctreeBranch>, // index -> same as normals but with branch_mask
    pub normal_mask: u8, // this space was being wasted anyway // if set, then this voxel has normal
    pub normals: Vec<Vec3d>, // index -> the index of the voxel in the normal_mask while ignoring 0's,
    // for eg, 01001000(normal_map), normal for 01000000(voxel) is 1 and for 00001000(voxel) is 0
    pub transparency_mask: u8,
    pub materials: Vec<u16>, // index from Octree.materials
}

impl VoxelOctreeBranch { // all branches consider their space as -1 to 1
    fn new(pos: OctreePos) -> Self {
        Self {
            pos, child_mask: 0, branch_mask: 0, chilranches: Vec::<Self>::new(), normal_mask: 0,
            normals: Vec::<Vec3d>::new(), materials: Vec::<u16>::new(), transparency_mask: 0,
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

    fn modify_voxel(&mut self, pos_mask: u8, material_index: u16, normal: Option<Vec3d>, transparency: bool) {
        let index = self.get_info_index(self.child_mask, pos_mask);
        self.materials[index] = material_index;

        if transparency && !(self.transparency_mask & pos_mask > 0) {
            self.transparency_mask = self.transparency_mask | pos_mask;
        } else if !transparency && self.transparency_mask & pos_mask > 0 {
            self.transparency_mask = self.transparency_mask & !pos_mask;
        }

        if let Some(normal) = normal { // if normal is supplied, either set it or change it
            if !(self.normal_mask & pos_mask > 0) {
                self.insert_normal(normal, pos_mask)
            } else {
                let index = self.get_info_index(self.normal_mask, pos_mask);
                self.normals[index] = normal;
            }

        } else if self.normal_mask & pos_mask > 0 { // no normal supplied + normal is there -> remove it
            let index = self.get_info_index(self.normal_mask, pos_mask);
            self.normals.remove(index);
            self.normal_mask = self.normal_mask & !pos_mask;
        }
    }

    fn set_voxel(&mut self, pos_mask: u8, material_index: u16, normal: Option<Vec3d>, transparency: bool) {
        if self.child_mask & pos_mask != 0 {
            self.modify_voxel(pos_mask, material_index, normal, transparency);
            return
        }
        self.child_mask = self.child_mask | pos_mask;
        
        if let Some(normal) =  normal {
            self.insert_normal(normal, pos_mask);
        }

        if transparency {
            self.transparency_mask = self.transparency_mask | pos_mask;
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

    fn insert_voxel_from_point(&mut self, point: Vec3d, depth: usize, material_index: u16, normal: Option<Vec3d>, transparency: bool) {
        let pos = self.get_pos_from_point(&point);
        let pos_mask = pos.get_mask();
        // println!("{:?}", pos);
        if depth == 0 {
            self.set_voxel(pos_mask, material_index, normal, transparency);
            return
        }
        let branch = self.add_or_get_branch(pos_mask, material_index, normal, transparency);
        branch.insert_voxel_from_point((point-branch.get_coord())*2.0, depth-1, material_index, normal, transparency);
    }
    
    fn add_or_get_branch(&mut self, pos_mask: u8, material_index: u16, normal: Option<Vec3d>, transparency: bool) -> &mut Self {
        if (self.branch_mask & pos_mask) > 0 {return self.get_branch_mut(pos_mask)}
        // if this isnt a branch
        self.add_branch(pos_mask, material_index, normal, transparency)
    }

    fn add_branch(&mut self, pos_mask: u8, material_index: u16, normal: Option<Vec3d>, transparency: bool) -> &mut Self {
        // create a new voxel
        if !(self.child_mask & pos_mask > 0) {
            self.set_voxel(pos_mask, material_index, normal, transparency);
        }

        // create a new branch
        self.branch_mask = self.branch_mask | pos_mask;
        let chilranch = Self::new(OctreePos::new(pos_mask));

        // insert branch at correct index
        let index = self.get_info_index(self.branch_mask, pos_mask);
        self.chilranches.insert(index, chilranch);
        &mut self.chilranches[index]
    }

    // find where the normal should go (index), and correctly set/modify normals and indices
    #[inline(always)]
    fn insert_normal(&mut self, normal: Vec3d, pos_mask: u8) {
        self.normal_mask = self.normal_mask | pos_mask;
        let index = self.get_info_index(self.normal_mask, pos_mask);
        self.normals.insert(index, normal);
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

    #[inline(always)]
    pub fn get_normal(&self, pos_mask: u8) -> &Vec3d {
        &self.normals[self.get_info_index(self.normal_mask, pos_mask)]
    }
    
    #[inline(always)]
    pub fn get_material_index(&self, child_mask: u8) -> u16 {
        self.materials[self.get_info_index(self.child_mask, child_mask)]
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

    #[inline]
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
    
    #[inline(always)]
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
    #[inline(always)]
    pub fn new(t: f64, plane: u8)-> Self {
        Self {t, plane}
    }

    // sort all ts in ascending order + remove any ts < current t
    #[inline(always)]
    pub fn get_next_hits(&self, mut hits: Vec<Self>) -> Vec<Self> {
        hits.sort_by(|a, b| a.t.partial_cmp(&b.t).unwrap());
        hits.iter().filter(|a| a.t > self.t).cloned().collect()
    }

    // same as above but only for the next emelent
    #[inline(always)]
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

    /// gives the normal of self.plane , direction is outside cube (dosent work when inside the cube (for now atleast))
    /// this depends on how the octree.hit handles planes
    #[inline(always)]
    pub fn get_plane_normal(&self) -> Vec3d {
        match self.plane {
            0b_10011001 => Vec3d::new(-1.0, 0.0, 0.0), // r -> normal towards left
            0b_01100110 => Vec3d::new(1.0, 0.0, 0.0), // l
            0b_00001111 => Vec3d::new(0.0, -1.0, 0.0), // u
            0b_11110000 => Vec3d::new(0.0, 1.0, 0.0), // d
            0b_11000011 => Vec3d::new(0.0, 0.0, -1.0), // b
            0b_00111100 => Vec3d::new(0.0, 0.0, 1.0), // f
            _ => panic!(),
        }
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
