
use rand::distributions::{Uniform, Distribution};
use std::f64::consts::PI;

use super::vec3d::Vec3d;

use super::world::{World, Camera};
use super::sphere::Sphere;
use super::material::{Material, Lambertian, Metal, Dielectric};


pub fn gen_world() -> World {

    let width = 720;
    let height = 480;
    let fov = PI/3.0;
    let samples_per_pixel = 100;
    let aperture = 0.1;
    let cam_position = Vec3d::new(0.0, 1.5, 0.0);
    let look_at = Vec3d::new(0.0, 1.0, -10.0);

    let mut world = World {
        cam: Camera::new(width, height, fov, samples_per_pixel, aperture, cam_position, look_at),
        objects: Vec::<Sphere>::new(),
    };

    world.objects.push( // ground
        Sphere {
            center: Vec3d::new(0.0, -1000.0, 0.0),
            radius: 1000.0,
            material: Material::Lambertian(
                Lambertian {
                    color:Vec3d::new(0.5, 0.5, 0.5)
                }
            ),
        }
    );

    world.objects.push(
        Sphere {
            center: Vec3d::new(0.0, 1.0, -10.0),
            radius: 1.0,
            material: Material::Metal(
                Metal {
                    color: Vec3d::new(0.77, 1.0, 0.77),
                    fuzz: 0.0,
                }
            ),
        }
    );
    world.objects.push(
        Sphere {
            center: Vec3d::new(-2.0, 1.0, -10.0),
            radius: 1.0,
            material: Material::Lambertian(
                Lambertian {
                    color: Vec3d::new(0.65, 1.0, 0.32),
                }
            ),
        }
    );
    world.objects.push(
        Sphere {
            center: Vec3d::new(2.0, 1.0, -10.0),
            radius: 1.0,
            material: Material::Dielectric(
                Dielectric {
                    // color: Vec3d::new(0.44, 0.21, 1.0),
                    color: Vec3d::new(1.0, 1.0, 1.0),
                    refractive_index: 1.5,
                    fuzz: 0.0,
                }
            ),
        }
    );
    world.objects.push( // trying hollow glass sphere
        Sphere {
            center: Vec3d::new(2.0, 1.0, -10.0),
            radius: -0.5,
            material: Material::Dielectric(
                Dielectric {
                    // color: Vec3d::new(0.44, 0.21, 1.0),
                    color: Vec3d::new(1.0, 1.0, 1.0),
                    refractive_index: 1.5,
                    fuzz: 0.0,
                }
            ),
        }
    );

    let mut rng = rand::thread_rng();
    let random = Uniform::new(0.0, 1.0);
    for _ in 0..70 {
        let mut z = (random.sample(&mut rng)-0.5)*20.0;
        z += -10.0;
        let x = (random.sample(&mut rng)-0.5)*20.0;
        let ground_radius = world.objects[0].radius;
        let ground_center = world.objects[0].center;
        let randomiser = random.sample(&mut rng);
        if randomiser < 0.333 {
            world.objects.push(
                Sphere {
                    center: (Vec3d::new(x, 0.0, z) - ground_center).unit()*(ground_radius+0.4) + ground_center,
                    radius: 0.4,
                    material: Material::Lambertian(
                        Lambertian {
                            color: Vec3d::new(random.sample(&mut rng), random.sample(&mut rng), random.sample(&mut rng)),
                        }
                    ),
                }
            );    
        } else if randomiser < 0.666 {
            world.objects.push(
                Sphere {
                    center: (Vec3d::new(x, 0.0, z) - ground_center).unit()*(ground_radius+0.4) + ground_center,
                    radius: 0.4,
                    material: Material::Metal(
                        Metal {
                            color: Vec3d::new(random.sample(&mut rng), random.sample(&mut rng), random.sample(&mut rng)),
                            fuzz: random.sample(&mut rng),
                        }
                    ),
                }
            );    
        } else {
            world.objects.push(
                Sphere {
                    center: (Vec3d::new(x, 0.0, z) - ground_center).unit()*(ground_radius+0.4) + ground_center,
                    radius: 0.4,
                    material: Material::Dielectric(
                        Dielectric {
                            color: Vec3d::new(random.sample(&mut rng), random.sample(&mut rng), random.sample(&mut rng)),
                            refractive_index: 1.0 + random.sample(&mut rng),
                            fuzz: random.sample(&mut rng)*0.0,
                        }
                    ),
                }
            );    
        }
    }

    world
}
