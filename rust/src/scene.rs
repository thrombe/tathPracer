
use rand::distributions::{Uniform, Distribution};
use std::f64::consts::PI;

use super::vec3d::Vec3d;

use super::world::{World, Camera};
use super::objects::{Object, Sphere, Plane, Triangle};
use super::material::{Material, Lambertian, Metal, Dielectric, Lit};


pub fn gen_world() -> World {

    let width = 720;
    let height = 480;
    let fov = PI/3.0;
    let samples_per_pixel = 100;
    let aperture = 0.1;
    let cam_position = Vec3d::new(-5.0, 2.50, 0.0);
    let look_at = Vec3d::new(0.0, 1.0, -10.0);

    let mut world = World {
        cam: Camera::new(width, height, fov, samples_per_pixel, aperture, cam_position, look_at),
        objects: Vec::<Object>::new(),
    };

    world.objects.push( // ground
        Object::Plane(Plane {
            point: Vec3d::new(0.0, 0.0, 0.0),
            normal: Vec3d::new(0.0, 1.0, 0.0),
            material: Material::Lambertian(Lambertian {
            // material: Material::Metal(Metal {
                color:Vec3d::new(0.45, 0.3, 0.45),
                // fuzz: 0.01,
            }),
        })
    );
    
    // world.objects.push( // sky
    //     Object::Sphere(Sphere {
    //         center: Vec3d::new(0.0, 0.0, 0.0),
    //         radius: 1000.0,
    //         material: Material::Lit(Lit {
    //             // color: Vec3d::new(0.9, 0.9, 9.5), // tinted sky
    //             color: Vec3d::new(0.0, 0.0, 0.0), // night
    //         }),
    //     })
    // );
    // world.objects.push( // sun(ish)
    //     Object::Sphere(Sphere {
    //         center: Vec3d::new(100.0, 100.0, 0.0),
    //         radius: 30.0,
    //         material: Material::Lit(Lit {
    //             color: Vec3d::new(50.0, 50.0, 50.0),
    //         }),
    //     })
    // );

    // world.objects.push( // tilted plane
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

    world.objects.push( // mid top
        Object::Sphere(Sphere {
            center: Vec3d::new(0.0, 3.0, -10.0),
            radius: 1.0,
            material: Material::Lit(Lit {
                color: Vec3d::new(2.80, 0.8, 3.2),
            }),
        })
    );
    world.objects.push( // mid
        Object::Sphere(Sphere {
            center: Vec3d::new(0.0, 1.0, -10.0),
            radius: 1.0,
            material: Material::Metal(Metal {
                color: Vec3d::new(0.77, 1.0, 0.77),
                fuzz: 0.0,
            }),
        })
    );
    world.objects.push( // left
        Object::Sphere(Sphere {
            center: Vec3d::new(-2.0, 1.0, -10.0),
            radius: 1.0,
            material: Material::Lambertian(Lambertian {
                color: Vec3d::new(0.65, 1.0, 0.32),
            }),
        })
    );
    world.objects.push( // left top
        Object::Sphere(Sphere {
            center: Vec3d::new(-2.0, 3.0, -10.0),
            radius: 1.0,
            material: Material::Metal(Metal {
                color: Vec3d::new(0.65, 0.23, 0.72),
                fuzz: 0.5,
            }),
        })
    );
    world.objects.push( // right
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
    world.objects.push( // trying hollow glass sphere
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
            world.objects.push(
                Object::Sphere(Sphere {
                    center: center,
                    radius: radius,
                    material: Material::Lambertian(Lambertian {
                        color: color,
                    }),
                })
            );    
        } else if randomiser < threshold*2.0 {
            world.objects.push(
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
            world.objects.push(
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
            world.objects.push(
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
