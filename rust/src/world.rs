

use super::vec3d::Vec3d;
use super::progress_indicator::ProgressIndicator;
use super::img;
// use super::math;

use rand::rngs::StdRng;
use rand::SeedableRng; // for rng
use rand::distributions::{Uniform, Distribution};

struct Camera {
    width: usize,
    height: usize,
    pos: Vec3d,
    fwd: Vec3d,
    up: Vec3d,
    right:Vec3d,
    fov: f64,
    scr_dist: f64,
    bouncy_depth: usize,
    t_correction: f64,
    far_away: f64,

    // focus_dist: f64,
    // aperture: f64,
}

fn gen_objects() -> Vec<Sphere> {
    let mut objects = Vec::<Sphere>::new();
    objects.push( // Ground
        Sphere {
            center: Vec3d::new(0.0, -1000.0, 0.0),
            radius: 1000.0,
            material: Material::Lambertian,
            color: Vec3d::new(0.5, 0.5, 0.5),
        }
    );

    objects.push(
        Sphere {
            center: Vec3d::new(0.0, 1.0, -10.0),
            radius: 1.0,
            material: Material::Metal,
            color: Vec3d::new(0.77, 1.0, 0.77),
        }
    );
    objects.push(
        Sphere {
            center: Vec3d::new(-2.0, 1.0, -10.0),
            radius: 1.0,
            material: Material::Lambertian,
            color: Vec3d::new(0.65, 1.0, 0.32),
        }
    );
    objects.push(
        Sphere {
            center: Vec3d::new(2.0, 1.0, -10.0),
            radius: 1.0,
            material: Material::Lambertian,
            color: Vec3d::new(0.44, 0.21, 1.0),
        }
    );

    let mut rng = rand::thread_rng();
    let random = Uniform::new(-1.0, 1.0);
    for _ in 0..70 {
        let z = -(random.sample(&mut rng)+1.0)*10.0;
        let x = (random.sample(&mut rng))*10.0;
        objects.push(
            Sphere {
                center: (Vec3d::new(x, 0.0, z) - objects[0].center).unit()*(objects[0].radius+0.4) + objects[0].center,
                radius: 0.4,
                material: Material::Lambertian,
                color: Vec3d::new(random.sample(&mut rng), random.sample(&mut rng), random.sample(&mut rng)),
            }
        );
    }

    objects
}

type Rng = rand::prelude::ThreadRng;

pub fn run_world() {
    let mut world = World {
        cam: Camera {
            width: 720,
            height: 480,
            fov: std::f64::consts::PI/3.0,
            pos: Vec3d::new(0.0, 1.5, 0.0),
            fwd: Vec3d::new(0.0, -0.0, -1.0),
            up: Vec3d::new(0.0, 1.0, 0.0),
            right: Vec3d::new(1.0, 0.0, 0.0),
            scr_dist: 0.0,
            bouncy_depth: 100,
            t_correction: 0.0000001,
            far_away: 1000000000.0,
        },
        objects: gen_objects(),
    };
    world.cam.scr_dist = {
        let cot = 1.0/(world.cam.fov/2.0).atan();
        ((world.cam.width as f64)/2.0)*cot
    };
    let samples = 10;
    let topleft = world.cam.fwd*world.cam.scr_dist 
                + world.cam.up*(world.cam.height as f64/2.0) 
                + world.cam.right*(-1.0)*(world.cam.width as f64/2.0);
    
    // rand stuff
    let half_pix_width = 1.0/(world.cam.width as f64*2.0);
    let half_pix_height = 1.0/(world.cam.height as f64*2.0);
    let rand_width_off = Uniform::new(-half_pix_width, half_pix_width);
    let rand_height_off = Uniform::new(-half_pix_height, half_pix_height);
    let mut rng = rand::thread_rng();
    // let mut rng = StdRng::from_entropy();

    let mut indicator = ProgressIndicator::new(world.cam.height);
    let mut img = img::new_img(world.cam.width as u32, world.cam.height as u32);
    for y in 0..world.cam.height as u32 {
        for x in 0..world.cam.width as u32 {
            let mut color = Vec3d::new(0.0, 0.0, 0.0);
            let ray = Ray::new(
                    world.cam.pos, 
                    topleft + world.cam.right*(x as f64) + world.cam.up*(-1.0)*(y as f64),
            );
            for s in 0..samples {
                let wiggle = world.cam.right*rand_width_off.sample(&mut rng) 
                           + world.cam.up*rand_height_off.sample(&mut rng);
                let mut ray = Ray::new(ray.pos, ray.dir + wiggle);
                color += world.get_ray_color(&mut ray, world.cam.bouncy_depth, &mut rng);
            }
            color *= 1.0/(samples as f64);
            color.x = color.x.sqrt();
            color.y = color.y.sqrt();
            color.z = color.z.sqrt();
            color *= 255.0;
            img::set(&mut img, x, y, color.x, color.y, color.z); // casting from f64 to u8 chops it at ends -> [0, 255]
        }
        indicator.indicate(y as usize);
    }
    img::dump_img(img);
}

struct Ray {
    pos: Vec3d,
    dir: Vec3d,
}

impl Ray {
    fn new(pos: Vec3d, dir: Vec3d) -> Self {
        Ray {pos, dir}
    }

    fn at(&mut self, t: f64) {
        self.pos = self.pos + self.dir*t;
    }
}

struct World {
    objects: Vec<Sphere>,
    cam: Camera,
}

impl World {
    fn hit(&self, ray: &Ray) -> (Option<&Sphere>, f64) {
        let mut min_t = self.cam.far_away;
        let mut hit_what = None;
        for object in &self.objects {
            if let Some(t) = object.hit(&ray, self.cam.t_correction) {
                if min_t > t {
                    min_t = t;
                    hit_what = Some(&*object);
                }
            }
        }
        (hit_what, min_t)
    }

    fn get_ray_color(&self, ray: &mut Ray, bouncy_depth: usize, rng: &mut Rng) -> Vec3d {
        if bouncy_depth <= 0 {return Vec3d::new(0.0, 0.0, 0.0)}

        if let (Some(hit_what), t) = self.hit(&ray) {
            ray.at(t);
            if let Some(mut ray) = hit_what.scatter(&ray, rng) {
                return self.get_ray_color(&mut ray, bouncy_depth-1, rng).mul(hit_what.color)
            }
            // ray absorbed
            return Vec3d::new(0.0, 0.0, 0.0)
        }

        // background color
        let t = (ray.dir.unit().y + 1.0)*0.5;
        return Vec3d::new(0.5, 0.7, 1.0).lerp(Vec3d::new(1.0, 1.0, 1.0), t)
    }
}

enum Material {
    Metal,
    Dielectric,
    Lambertian,
}

impl Material {
    fn scatter(&self, ray: &Ray, normal: &Vec3d, rng: &mut Rng) -> Option<Ray> {
        let rand = Uniform::new(-1.0, 1.0);
        match self {
            Material::Lambertian => {
                let mut dir = *normal + Vec3d::new(rand.sample(rng), rand.sample(rng), rand.sample(rng));
                // let tea = 0.00000001;
                // if dir.x < tea && dir.y < tea && dir.z < tea { // ray dir goes 0
                //     dir = *normal;
                // }

                Some(Ray::new(ray.pos, dir))
            },
            Material::Metal => {
                let a = ray.dir.unit();
                let reflected = a - *normal*2.0*a.dot(*normal);
                if normal.dot(reflected) < 0.0 {return None}
                Some(Ray::new(ray.pos, reflected))
            },
            Material::Dielectric => {
                None
            },
        }
    }
}

struct Sphere {
    center: Vec3d,
    radius: f64,
    material: Material,
    // instead of material being u8, i can try using enum or something for this
    // or, matrial points to another struct, that receivers an immutable pointer
    // to sphere, and thus can receive normal and stuff. the material struct just
    // needs to give next scattered ray and give some color to the ray.
    color: Vec3d, // maybe store this in material?
}

pub fn test() {
    let sp = Sphere {
        center: Vec3d::new(0.0, 0.0, -5.0),
        radius: 1.0,
        material: Material::Lambertian,
        color: Vec3d::new(1.0, 0.0, 0.0),
    };
    let ray = Ray {
        pos: Vec3d::new(0.0, 0.0, 0.0),
        dir: Vec3d::new(0.0, 0.0, -1.0),
    };
    let val = sp.hit(&ray, 0.0001);
    println!("{:?}", val);
}

impl Sphere {
    /// returns t for the ray if hit, else returns None
    fn hit(&self, ray: &Ray, t_correction: f64) -> Option<f64> {
        // if hit return Some(t)
        // otherwise return None
        let oc = ray.pos-self.center;
        let neg_b = -ray.dir.dot(oc);
        let b_sq = ray.dir.dot(ray.dir);
        let d_by_4 = neg_b*neg_b - b_sq*(oc.dot(oc) - self.radius*self.radius);
        if d_by_4 < 0.0 { // didnt hit
            return None
        }
        let t = (neg_b - d_by_4.sqrt())/b_sq;
        if t < t_correction { // ray hitting behind the camera or really close to object ( t < 0.0000001 for really close thing)
            return None
        }
        Some(t)
    }

    fn scatter(&self, ray: &Ray, rng: &mut Rng) -> Option<Ray> {
        self.material.scatter(ray, &self.normal(&ray.pos), rng)
    }

    fn normal(&self, point: &Vec3d) -> Vec3d {
        (*point - self.center).unit()
    }
}
