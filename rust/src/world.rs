

use super::vec3d::Vec3d;
use super::progress_indicator::ProgressIndicator;
use super::img;

use super::ray::Ray;
use super::scene::gen_world;
use super::sphere::Sphere;

// use rand::rngs::StdRng;
// use rand::SeedableRng; // for rng
use rand::distributions::{Uniform, Distribution};

pub struct Camera {
    pub width: usize,
    pub height: usize,
    pos: Vec3d,
    fwd: Vec3d,
    up: Vec3d,
    right:Vec3d,
    fov: f64,
    pub bouncy_depth: usize,
    pub samples_per_pixel: usize,
    pub t_correction: f64,
    pub far_away: f64,
    topleft: Vec3d,
    scr_dist: f64,
}

impl Camera {

    pub fn new(width: usize, height: usize, fov: f64, samples_per_pixel: usize) -> Camera {
        let mut cam = Camera {
            // inputs
            width, height, fov, samples_per_pixel,
            
            // calculated
            topleft: Vec3d::zero(), scr_dist: 0.0,
            
            // defaults
            pos: Vec3d::new(0.0, 0.0, 0.0),
            fwd: Vec3d::new(0.0, 0.0, -1.0),
            up: Vec3d::new(0.0, 1.0, 0.0),
            right: Vec3d::new(1.0, 0.0, 0.0),
            t_correction: 0.0000001,
            far_away: 1000000000.0,
            bouncy_depth: 100,
        };
        cam.setup();
        cam
    }

    pub fn setup(&mut self) {
        self.scr_dist = {
                let cot = 1.0/(self.fov/2.0).atan();
                ((self.width as f64)/2.0)*cot
            };
        self.topleft = self.pos 
                     + self.fwd*self.scr_dist
                     + self.up*(self.height as f64/2.0) 
                     + self.right*(-1.0)*(self.width as f64/2.0);
    }

    pub fn get_ray(&self, x: u32, y: u32) -> Ray {
        Ray::new(self.pos, self.topleft + self.right*(x as f64) + self.up*(-1.0)*(y as f64))
    }

    pub fn move_to(&mut self, pos: Vec3d) {
        self.pos = pos;

        self.setup();
    }

    /// this works as long as the camera is kind of vertical (i.e. not upside down or something)
    pub fn face_at(&mut self, dir: Vec3d) {
        self.fwd = dir.unit();
        self.right = self.fwd.cross(self.up);
        self.up = self.right.cross(self.fwd);

        self.setup();
    }
}

pub type Rng = rand::prelude::ThreadRng;

pub fn run_world() {
    let world = gen_world();
    let samples = world.cam.samples_per_pixel;
    
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
            for _ in 0..samples {
                let mut ray = world.cam.get_ray(x, y);
                let wiggle = world.cam.right*rand_width_off.sample(&mut rng) 
                           + world.cam.up*rand_height_off.sample(&mut rng);
                ray.dir += wiggle;
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


pub struct World {
    pub objects: Vec<Sphere>,
    pub cam: Camera,
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
            ray.new_pos(t);
            if let Some(mut ray) = hit_what.scatter(&ray, rng) {
                return self.get_ray_color(&mut ray, bouncy_depth-1, rng).mul(*hit_what.color())
            }
            // ray absorbed
            return Vec3d::new(0.0, 0.0, 0.0)
        }

        // background color
        let t = (ray.dir.unit().y + 1.0)*0.5;
        return Vec3d::new(0.5, 0.7, 1.0).lerp(Vec3d::new(1.0, 1.0, 1.0), t)
    }
}

