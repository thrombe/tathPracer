

use super::vec3d::Vec3d;
use super::progress_indicator::ProgressIndicator;
use super::img;

use super::ray::{Ray, RayHitfo};
use super::scene::gen_world;
use super::objects::{Object};
use super::material::Material;

// use rand::rngs::StdRng;
// use rand::SeedableRng; // for rng
use rand::distributions::{Uniform, Distribution};
use std::thread; // for multi-threading
use std::sync::Arc;

pub struct Camera {
    // public for access from world
    // parameters
    pub width: usize,
    pub height: usize,
    pub bouncy_depth: usize,
    pub samples_per_pixel: usize,
    pub t_correction: f64,
    pub far_away: f64,
    fov: f64,
    lens_radius: f64,
    pos: Vec3d,
    look_at: Vec3d,

    // cache
    focus_length: f64,
    topleft: Vec3d,
    scr_dist: f64,
    scr_horizontal: Vec3d,
    scr_vertical: Vec3d,

    fwd: Vec3d,
    up: Vec3d,
    right:Vec3d,
    rng_distribution: Uniform<f64>,
    cores: usize,
}

impl Camera {

    pub fn new(width: usize, height: usize, fov: f64, samples_per_pixel: usize, aperture: f64, pos: Vec3d, look_at: Vec3d) -> Camera {
        let mut cam = Camera {
            // inputs
            width, height, fov, samples_per_pixel, lens_radius: aperture/2.0, pos, look_at,
            
            // calculated
            topleft: Vec3d::zero(), scr_horizontal: Vec3d::zero(), scr_vertical: Vec3d::zero(), scr_dist: 0.0, focus_length: 0.0,
            
            // defaults that should not be touched
            fwd: Vec3d::new(0.0, 0.0, -1.0),
            up: Vec3d::new(0.0, 1.0, 0.0),
            right: Vec3d::new(1.0, 0.0, 0.0),
            rng_distribution: Uniform::new(-1.0, 1.0),
            
            // defaults
            t_correction: 0.0000001,
            far_away: 1000000000.0,
            bouncy_depth: 100,
            cores: 7,
        };
        cam.setup();
        cam
    }

    pub fn setup(&mut self) {
        self.fwd = self.look_at-self.pos;
        self.focus_length = self.fwd.size();
        self.fwd *= 1.0/self.focus_length;
        self.right = self.fwd.cross(self.up);
        self.up = self.right.cross(self.fwd);


        self.scr_dist = {
                let cot = 1.0/(self.fov/2.0).tan();
                ((self.width as f64)/2.0)*cot
            };
        
        self.scr_horizontal = self.right*(self.focus_length/self.scr_dist);
        self.scr_vertical = self.up*(self.focus_length/self.scr_dist)*(-1.0); // -1 to reduce calc in get_ray
        self.topleft = self.fwd*self.focus_length // topleft is direction only(relative to cam), not exact coords
                     + self.scr_vertical*(-1.0)*(self.height as f64/2.0)
                     + self.scr_horizontal*(-1.0)*(self.width as f64/2.0);
    }

    pub fn get_ray(&self, x: f64, y: f64, rng: &mut Rng) -> Ray {
        let random_in_unit_disc = {
            let mut p: (f64, f64);
            loop {
                p = (self.rng_distribution.sample(rng), self.rng_distribution.sample(rng));
                if p.0*p.0+p.1*p.1 > 1.0 {continue}
                break
            }
            p
        };
        let mut offset = self.right*random_in_unit_disc.0 + self.up*random_in_unit_disc.1;
        offset *= self.lens_radius;
        // this ray formation is like shifting the camera by some offset, while looking at the same point
        Ray::new(self.pos + offset, self.topleft - offset + self.scr_horizontal*x + self.scr_vertical*y)
    }

    pub fn set_focus_length(&mut self, focus_length: f64) {
        self.focus_length = focus_length;

        self.setup();
    }
}

pub type Rng = rand::prelude::ThreadRng;

pub fn run_world() {
    let world = Arc::new(gen_world()); // there should be no arc penalty for use in hot loops as the penalty of rc is only when refrences get created/destroyed. right??
    let samples = world.cam.samples_per_pixel;
    let cores = world.cam.cores;
    
    let samples_per_thread = samples/cores;
    let cores_without_extra_work = {
        let leftover_work = samples-samples_per_thread*cores;
        cores-leftover_work
    };

    let mut worker_vec = vec![];
    let mut img_buffer: Vec<Vec3d> = Vec::new();
    for core in 0..cores {
        let world = Arc::clone(&world);
        let samples = if core >= cores-cores_without_extra_work {samples_per_thread} else {samples_per_thread+1};
        let mut indicator = if core == cores-1 {ProgressIndicator::new(world.cam.height)} else {ProgressIndicator::new(0)};

        let mut process = move || -> Vec<Vec3d> { // cant move this outside the loop as rust complains that "world" dies earlier than this closure
            // rand stuff
            let rand_off = Uniform::new(-0.5, 0.5);
            let mut rng = rand::thread_rng();
                
            let mut img_buffer: Vec<Vec3d> = Vec::with_capacity(world.cam.width*world.cam.height);
    
            for y in 0..world.cam.height as u32 {
                for x in 0..world.cam.width as u32 {
                    let mut color = Vec3d::new(0.0, 0.0, 0.0);
                    for _ in 0..samples {
                        let mut ray = world.cam.get_ray(
                            x as f64 + rand_off.sample(&mut rng), 
                            y as f64 + rand_off.sample(&mut rng), 
                            &mut rng,
                            );
                        color += world.get_ray_color(&mut ray, world.cam.bouncy_depth, &mut rng);
                    }
                    img_buffer.push(color);
                }
                indicator.indicate(y as usize);
            }
            img_buffer
        };    

        if core == cores-1 {
            img_buffer = process();
        } else {
            let worker = thread::spawn(process);
            worker_vec.push(worker);
        }
    }
    for worker in worker_vec {
        if let Ok(buffer) = worker.join() {
            for i in 0..buffer.len() {
                img_buffer[i] += buffer[i];
            }
        }
    }

    let mut img = img::new_img(world.cam.width as u32, world.cam.height as u32);
    for y in 0..world.cam.height {
        for x in 0..world.cam.width {
            let buffer_index = y*world.cam.width + x;
            let mut color = img_buffer[buffer_index];
            color *= 1.0/(world.cam.samples_per_pixel as f64);
            color.x = color.x.sqrt();
            color.y = color.y.sqrt(); // cant do this in parallel as sqrt(x+y) != sqrt(x)+sqrt(y)
            color.z = color.z.sqrt();
            color *= 255.0;
            img::set(&mut img, x as u32, y as u32, color.x, color.y, color.z); // casting from f64 to u8 chops it at ends -> [0, 255]
        }
    }

    img::dump_img(img);
}


pub struct World {
    pub objects: Vec<Object>, // Vec<Box<dyn hittable>> or enum ?
    pub cam: Camera,
}

impl World {
    #[inline(always)]
    fn hit(&self, ray: &Ray) -> Option<RayHitfo> {
        let mut min_t = self.cam.far_away;
        let mut hit_what = None;
        for object in &self.objects {
            if let Some(hitfo) = object.hit(&ray, self.cam.t_correction) {
                if min_t > hitfo.t {
                    min_t = hitfo.t;
                    hit_what = Some(hitfo);
                }
            }
        }
        hit_what
    }

    fn get_ray_color(&self, ray: &mut Ray, bouncy_depth: usize, rng: &mut Rng) -> Vec3d {
        if bouncy_depth <= 0 {return Vec3d::new(0.0, 0.0, 0.0)}

        if let Some(hitfo) = self.hit(&ray) {
            // lit objects
            if let Material::Lit(obj) = hitfo.material {
                return obj.color
            }

            if let Some(mut ray) = hitfo.scatter(rng) {
                return self.get_ray_color(&mut ray, bouncy_depth-1, rng).mul(*hitfo.material.color())
            }
            // ray absorbed
            return Vec3d::new(0.0, 0.0, 0.0)
        }

        // background color
        let t = (ray.dir.unit().y + 1.0)*0.5;
        return Vec3d::new(0.5, 0.7, 1.0).lerp(Vec3d::new(1.0, 1.0, 1.0), t)
    }
}

