

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
    pub pos: Vec3d,
    pub fwd: Vec3d,
    pub up: Vec3d,
    pub right:Vec3d,
    pub fov: f64,
    pub scr_dist: f64,
    pub bouncy_depth: usize,
    pub samples_per_pixel: usize,
    pub t_correction: f64,
    pub far_away: f64,

    // focus_dist: f64,
    // aperture: f64,
}

pub type Rng = rand::prelude::ThreadRng;

pub fn run_world() {
    let world = gen_world();
    let samples = world.cam.samples_per_pixel;
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
            for _ in 0..samples {
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

