
use std::fs::File;
use std::io::Read;

use super::triangle_octree::{Triangle, TriangleOctree};
use super::vec3d::Vec3d;

impl TriangleOctree {
    pub fn import_from_obj(&mut self, path: &str, material_index: u16) {
        let mut file = File::open(path).unwrap();
        let mut data = String::new();
        file.read_to_string(&mut data).unwrap();

        let mut triangles = vec!();
        for line in data.lines() {
            let line = line.split_whitespace().collect::<Vec<&str>>();
            if line.is_empty() {continue}
            match line[0] {
                "v" => self.vertices.push(Vec3d::new(
                    line[1].parse().unwrap(),
                    line[2].parse().unwrap(),
                    line[3].parse().unwrap())),
                "f" => {
                    let mut triangle = Triangle {
                        vertex_indices: (
                            line[1].parse::<u32>().unwrap()-1,
                            line[2].parse::<u32>().unwrap()-1,
                            line[3].parse::<u32>().unwrap()-1,    
                        ),
                        material_index,
                        normal: Vec3d::zero(),
                    };
                    let vertices = triangle.get_vertices(&self.vertices);
                    // triangle.normal = (*vertices.2-*vertices.0).cross(*vertices.1-*vertices.0).unit(); // clock
                    triangle.normal = (*vertices.1-*vertices.0).cross(*vertices.2-*vertices.0).unit(); // anti clock
                    triangles.push(triangle);
                }
                _ => continue
            }
        }

        for triangle in triangles {
            self.insert_triangle(triangle);
        }
    }
}