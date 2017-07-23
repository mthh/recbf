extern crate image;
extern crate recbf;

use std::path::Path;
use recbf::recursive_bf;

use image::GenericImage;

fn main() {
    // Open the image to use:
    let img = image::open(&Path::new("examples/in.jpg")).unwrap();
    // How many channel ?
    let channel: u32 = match img.color() {
        image::ColorType::Gray(_) => 1,
        image::ColorType::GrayA(_) => 2,
        image::ColorType::RGB(_) => 3,
        image::ColorType::RGBA(_) => 4,
        _ => panic!("I don't know how many channels!")
    };
    let (width, height) : (u32, u32) = img.dimensions();
    // Read the image content in a mutable variable:
    let mut content = img.raw_pixels();
    // The content will be modified inplace:
    recursive_bf(&mut content, 0.03, 0.1, width, height, channel);
    // Save the result:
    image::save_buffer(Path::new("examples/out.jpg"), &content, width, height, img.color());
}
