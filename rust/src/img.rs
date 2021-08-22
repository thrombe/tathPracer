
use image::{ImageBuffer, RgbImage};
use std::path::Path;
use std::time::{SystemTime, UNIX_EPOCH};

/// an example of how a simple image generation looks like
pub fn example_img() {
    let width: u32 = 255;
    let height: u32 = 255;
    let mut img = new_img(width, height);
    for y in 0..width {
        for x in 0..height {
            // set(&mut img, x, y, 255.0-x as f64, y as f64, x as f64);
            set_u8(&mut img, x, y, 255-x as u8, y as u8, x as u8);
        }
    }
    dump_img(img)
}

pub fn new_img(width: u32, height: u32) -> ImageBuffer<image::Rgb<u8>, std::vec::Vec<u8>> {
    let image: RgbImage = ImageBuffer::new(width, height);
    image
}

/// finds what file name is valid
fn file_name() -> String {
    let now: u64 = SystemTime::now()
        .duration_since(UNIX_EPOCH).unwrap()
        .as_secs();
    let path = format!("../images/{}.png", now);
    if Path::new(&path).exists() { // should never happen, but i dont wanna loose files
        println!("file already exists: {}", &path);
        return file_name()
    }
    path
}

/// dumps image with a valid filename
pub fn dump_img(img: ImageBuffer<image::Rgb<u8>, std::vec::Vec<u8>>) {
    img.save(file_name()).unwrap();
}

pub fn dump_img_mut(img: &mut ImageBuffer<image::Rgb<u8>, std::vec::Vec<u8>>) {
    img.save(file_name()).unwrap();
}

/// set pixels of image
#[inline(always)]
pub fn set(img: &mut ImageBuffer<image::Rgb<u8>, std::vec::Vec<u8>>, x: u32, y: u32, r: f64, g: f64, b: f64) {
    img.put_pixel(x, y, image::Rgb([r as u8, g as u8, b as u8]))
}

/// set pixels of image
#[inline(always)]
pub fn set_u8(img: &mut ImageBuffer<image::Rgb<u8>, std::vec::Vec<u8>>, x: u32, y: u32, r: u8, g: u8, b: u8) {
    img.put_pixel(x, y, image::Rgb([r, g, b]))
}

#[allow(non_camel_case_types)]
pub struct pix {
    pub x: u32,
    pub y: u32,
    pub r: f64,
    pub g: f64,
    pub b: f64,
}

impl pix {

    #[inline(always)]
    pub fn new(x: u32, y: u32, r: f64, g: f64, b: f64) -> pix {
        pix { x, y, r, g, b}
    }
}