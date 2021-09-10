
use std::fs::File;
use std::io::Read;

use super::triangle_octree::{Triangle, TriangleOctree};
use super::vec3d::Vec3d;
use super::math;

impl TriangleOctree {
    pub fn import_from_obj(&mut self, path: &str, material_index: u16) {
        let mut file = File::open(path).unwrap();
        let mut data = String::new();
        file.read_to_string(&mut data).unwrap();

        let mut triangles = vec!();
        for line in data.lines() {
            let mut line = line.split_whitespace().collect::<Vec<&str>>();
            if line.is_empty() {continue}
            match line[0] {
                "v" => self.vertices.push(Vec3d::new(
                    line[1].parse().unwrap(),
                    line[2].parse().unwrap(),
                    line[3].parse().unwrap())),
                "f" => {
                    line[1] = line[1].split("/").collect::<Vec<&str>>()[0];
                    line[2] = line[2].split("/").collect::<Vec<&str>>()[0];
                    line[3] = line[3].split("/").collect::<Vec<&str>>()[0];
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

        let mut max_len = 0.0;
        for vertex in &self.vertices {
            max_len = math::max(max_len, vertex.size());
        }
        max_len = 1.0/max_len;
        for i in 0..self.vertices.len() {
            self.vertices[i] = self.vertices[i]*max_len*self.half_size;
        }

        for triangle in triangles {
            self.insert_triangle(triangle);
        }
    }
}