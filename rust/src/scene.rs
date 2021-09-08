#![allow(unused_variables)]
#![allow(unused_imports)]

use rand::distributions::{Uniform, Distribution};
use std::f64::consts::PI;

use super::vec3d::Vec3d;

use super::world::{World, Camera};
use super::objects::{Object, Sphere, Plane};
use super::voxel_octree::{VoxelOctree};
use super::object_octree::ObjectOctree;
use super::triangle_octree::{Triangle, TriangleOctree};
use super::material::{Material, Lambertian, Metal, Dielectric, Lit};


pub fn gen_world() -> World {
    // spheres()
    // voxel_octree()
    triangle_octree()
}

fn spheres() -> World {

    let width = 720;
    let height = 480;
    let fov = PI/3.0;
    let samples_per_pixel = 100;
    let aperture = 0.1;
    let cam_position = Vec3d::new(-5.0, 2.50, 0.0);
    let look_at = Vec3d::new(0.0, 1.0, -10.0);

    let mut world = World {
        cam: Camera::new(width, height, fov, samples_per_pixel, aperture, cam_position, look_at),
        octree: ObjectOctree::new(Vec3d::zero(), 50.0),
    };
    
    world.octree.insert_object( // ground
        Object::Plane(Plane {
            point: Vec3d::new(0.0, 0.0, 0.0),
            normal: Vec3d::new(0.0, 1.0, 0.0),
            material: Material::Lambertian(Lambertian {
            // material: Material::Metal(Metal {
                color:Vec3d::new(0.45, 0.3, 0.45),
                // fuzz: 0.01
            }),
        })
    );
    
    // world.octree.insert_object( // sky
    //     Object::Sphere(Sphere {
    //         center: Vec3d::new(0.0, 0.0, 0.0),
    //         radius: 1000.0,
    //         material: Material::Lit(Lit {
    //             // color: Vec3d::new(0.9, 0.9, 9.5), // tinted sky
    //             color: Vec3d::new(0.0, 0.0, 0.0), // night
    //         }),
    //     })
    // );
    // world.octree.insert_object( // sun(ish)
    //     Object::Sphere(Sphere {
    //         center: Vec3d::new(100.0, 100.0, 0.0),
    //         radius: 30.0,
    //         material: Material::Lit(Lit {
    //             color: Vec3d::new(50.0, 50.0, 50.0),
    //         }),
    //     })
    // );

    // world.octree.insert_object( // tilted plane
    //     Object::Plane(Plane {
    //         point: Vec3d::new(-15.0, 0.0, 0.0),
    //         normal: Vec3d::new(1.0, -0.4, 0.6),
    //         // material: Material::Lambertian(Lambertian {
    //         //     color: Vec3d::new(0.7, 0.5, 0.5),
    //         // }),
    //         material: Material::Metal(Metal {
    //             color: Vec3d::new(1.0, 1.0, 1.0),
    //             fuzz: 0.0,
    //         })
    //     })
    // );

    world.octree.insert_object( // mid top
        Object::Sphere(Sphere {
            center: Vec3d::new(0.0, 3.0, -10.0),
            radius: 1.0,
            material: Material::Lit(Lit {
                color: Vec3d::new(2.80, 0.8, 3.2),
            }),
        })
    );
    world.octree.insert_object( // mid
        Object::Sphere(Sphere {
            center: Vec3d::new(0.0, 1.0, -10.0),
            radius: 1.0,
            material: Material::Metal(Metal {
                color: Vec3d::new(0.77, 1.0, 0.77),
                fuzz: 0.0,
            }),
        })
    );
    world.octree.insert_object( // left
        Object::Sphere(Sphere {
            center: Vec3d::new(-2.0, 1.0, -10.0),
            radius: 1.0,
            material: Material::Lambertian(Lambertian {
                color: Vec3d::new(0.65, 1.0, 0.32),
            }),
        })
    );
    world.octree.insert_object( // left top
        Object::Sphere(Sphere {
            center: Vec3d::new(-2.0, 3.0, -10.0),
            radius: 1.0,
            material: Material::Metal(Metal {
                color: Vec3d::new(0.65, 0.23, 0.72),
                fuzz: 0.5,
            }),
        })
    );
    world.octree.insert_object( // right
        Object::Sphere(Sphere {
            center: Vec3d::new(2.0, 1.0, -10.0),
            radius: 1.0,
            material: Material::Dielectric(Dielectric {
                // color: Vec3d::new(0.44, 0.21, 1.0),
                color: Vec3d::new(1.0, 1.0, 1.0),
                refractive_index: 1.5,
                fuzz: 0.0,
            }),
        })
    );
    world.octree.insert_object( // trying hollow glass sphere
        Object::Sphere(Sphere {
            center: Vec3d::new(2.0, 1.0, -10.0),
            radius: -0.3,
            material: Material::Dielectric(Dielectric {
                // color: Vec3d::new(0.44, 0.21, 1.0),
                color: Vec3d::new(1.0, 1.0, 1.0),
                refractive_index: 1.5,
                fuzz: 0.0,
            }),
        })
    );

    let mut rng = rand::thread_rng();
    let random = Uniform::new(0.0, 1.0);
    for _ in 0..90 {
        let z = (random.sample(&mut rng)-0.5)*20.0 - 10.0;
        let x = (random.sample(&mut rng)-0.5)*20.0;
        let radius = 0.4;
        let center = Vec3d::new(x, radius, z);
        let randomiser = random.sample(&mut rng);
        let color = Vec3d::new(random.sample(&mut rng), random.sample(&mut rng), random.sample(&mut rng));
        let threshold = 1.0/4.0; // 4 types of spheres
        if randomiser < threshold {
            world.octree.insert_object(
                Object::Sphere(Sphere {
                    center: center,
                    radius: radius,
                    material: Material::Lambertian(Lambertian {
                        color: color,
                    }),
                })
            );    
        } else if randomiser < threshold*2.0 {
            world.octree.insert_object(
                Object::Sphere(Sphere {
                    center: center,
                    radius: radius,
                    material: Material::Metal(Metal {
                        color: color,
                        fuzz: random.sample(&mut rng),
                    }),
                })
            );    
        } else if randomiser < threshold*3.0 {
            world.octree.insert_object(
                Object::Sphere(Sphere {
                    center: center,
                    radius: radius,
                    material: Material::Dielectric(Dielectric {
                        color: color,
                        refractive_index: 1.0 + random.sample(&mut rng),
                        fuzz: random.sample(&mut rng)*0.0,
                    }),
                })
            );
        } else {
            world.octree.insert_object(
                Object::Sphere(Sphere {
                    center: center,
                    radius: radius,
                    material: Material::Lit(Lit {
                        color: color*4.0,
                    }),
                })
            );
        }
    }

    world
}

fn voxel_octree() -> World {

    let width = 720;
    let height = 480;
    let fov = PI/3.0;
    let samples_per_pixel = 100;
    let aperture = 0.0;
    let cam_position = Vec3d::new(1.3, 0.7, 4.0);
    let look_at = Vec3d::new(0.0, 0.0, 0.0);

    let mut world = World {
        cam: Camera::new(width, height, fov, samples_per_pixel, aperture, cam_position, look_at),
        octree: ObjectOctree::new(Vec3d::zero(), 100.0),
    };
    // world.cam.bouncy_depth = 1000;

    let mut oct = VoxelOctree::new(Vec3d::zero(), 2.0, Some(1));
    
    let material = Material::Lambertian(Lambertian {
        color: Vec3d::new(0.7, 0.4, 0.7),
    });
    let material_index = oct.add_material(material);

    // let material = Material::Lambertian(Lambertian {
    let material = Material::Dielectric(Dielectric {
        color: Vec3d::new(0.99, 0.9, 0.9),
        refractive_index: 1.3,
        fuzz: 0.0,
    });
    let material_index2 = oct.add_material(material);

    let mut rng = rand::thread_rng();
    let random = Uniform::new(-1.0, 1.0);
    // for _ in 0..100 {
    //     let point = Vec3d::new(random.sample(&mut rng), random.sample(&mut rng), 0.01);
    //     oct.insert_voxel(point, 1, material_index, None);
    // }
    for _ in 0..500 {
        let point = Vec3d::new(random.sample(&mut rng), random.sample(&mut rng), random.sample(&mut rng));
        let normal = point.clone().unit();
        oct.insert_voxel(point, 1, material_index2, None);
    }

    let point = Vec3d::new(-0.1, 0.01, -0.0);
    let normal = point.clone();
    oct.insert_voxel(point, 1, material_index, None);
    // let point = Vec3d::new(0.99, 0.99, 0.99);
    // let normal = point.clone();
    // oct.insert_voxel(point, 0, material_index2, None);

    world.octree.insert_object(Object::VoxelOctree(oct));

    // debug rays
    // use super::ray::Ray;
    // let mut rng = rand::thread_rng();
    // let ray = Ray::new(world.cam.pos, world.cam.fwd);
    // if let Some(hitfo) = world.objects[0].hit(&ray, world.cam.t_correction, &mut rng) {
    //     dbg!(&hitfo);
    //     if let Some(ray) = hitfo.material.scatter(&hitfo.ray, &hitfo.normal, &mut rng) {
    //         if let Some(hitfo) = world.objects[0].hit(&ray, world.cam.t_correction, &mut rng) {
    //             dbg!(&hitfo);
    //         }
    //     }
    // }

    world
}

fn triangle_octree() -> World {

    let width = 720;
    let height = 480;
    let fov = PI/3.0;
    let samples_per_pixel = 100;
    let aperture = 0.0;
    let cam_position = Vec3d::new(1.3, 0.7, 4.0);
    let look_at = Vec3d::new(0.0, 0.0, 0.0);

    let mut world = World {
        cam: Camera::new(width, height, fov, samples_per_pixel, aperture, cam_position, look_at),
        octree: ObjectOctree::new(Vec3d::zero(), 100.0),
    };

    let mut oct = TriangleOctree::new(Vec3d::zero(), 50.0);
    
    let material = Material::Lambertian(Lambertian {
        color: Vec3d::new(0.7, 0.4, 0.7),
    });
    let material_index = oct.add_material(material);

    oct.vertices.push(Vec3d::new(-1.0, 1.0, 0.0));
    oct.vertices.push(Vec3d::new(1.0, 1.0, 0.0));
    oct.vertices.push(Vec3d::new(0.0, -1.0, 0.0));
    oct.insert_triangle(Triangle {
        vertex_indices: (0, 1, 2),
        material_index: material_index,
        normal: Vec3d::new(0.0, 0.0, 1.0),
    });

    world.octree.insert_object(Object::TriangleOctree(oct));

    world
}
