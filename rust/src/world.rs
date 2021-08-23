

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
}

impl Camera {

    pub fn new(width: usize, height: usize, fov: f64, samples_per_pixel: usize, aperture: f64, pos: Vec3d, look_at: Vec3d) -> Camera {
        let mut cam = Camera {
            // inputs
            width, height, fov, samples_per_pixel, lens_radius: aperture/2.0, pos, look_at,
            
            // calculated
            topleft: Vec3d::zero(), scr_horizontal: Vec3d::zero(), scr_vertical: Vec3d::zero(), scr_dist: 0.0, focus_length: 0.0,
            
            // defaults
            fwd: Vec3d::new(0.0, 0.0, -1.0),
            up: Vec3d::new(0.0, 1.0, 0.0),
            right: Vec3d::new(1.0, 0.0, 0.0),
            t_correction: 0.0000001,
            far_away: 1000000000.0,
            bouncy_depth: 100,
            rng_distribution: Uniform::new(-1.0, 1.0),
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
    let world = gen_world();
    let samples = world.cam.samples_per_pixel;
    
    // rand stuff
    let rand_off = Uniform::new(-0.5, 0.5);
    let mut rng = rand::thread_rng();
    // let mut rng = StdRng::from_entropy();

    let mut indicator = ProgressIndicator::new(world.cam.height);
    let mut img = img::new_img(world.cam.width as u32, world.cam.height as u32);
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

