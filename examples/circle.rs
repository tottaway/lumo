use glam::DVec3;
use rust_tracer::tracer::*;
use rust_tracer::*;
use std::f64::consts::PI;

// s = v = 1.0. h in radians
fn hsv_to_rgb(h: f64) -> DVec3 {
    let f = |n: f64| {
        let k = (n + h / (PI / 3.0)) % 6.0;
        1.0 - k.min(4.0 - k).min(1.0).max(0.0)
    };

    let rgb = DVec3::new(f(5.0), f(3.0), f(1.0)) * 255.0;

    srgb_to_linear(rgb.x as u8, rgb.y as u8, rgb.z as u8)
}

fn main() -> Result<(), std::io::Error> {
    let camera = Camera::new(
        DVec3::new(0.0, 1.5, 1.5),
        DVec3::ZERO,
        DVec3::new(0.0, 1.0, -1.0),
        1.0,
    );

    let mut scene = Scene::default();
    let ground = -0.2;

    scene.add(Plane::new(
        DVec3::new(0.0, ground, 0.0),
        DVec3::new(0.0, 1.0, 0.0),
        Material::metal(Texture::Solid(srgb_to_linear(166, 44, 43)), 0.01),
    ));

    let r = 0.2;
    scene.add(Sphere::new(
        DVec3::new(0.0, ground + r + 0.1, 0.0),
        r,
        Material::Light(Texture::Solid(srgb_to_linear(255, 255, 255))),
    ));

    scene.add(Sphere::new(
        DVec3::ZERO,
        3.0,
        Material::diffuse(Texture::Solid(srgb_to_linear(0, 0, 0))),
    ));

    let circle_s = 8;
    let offset = PI / circle_s as f64;
    let r = 0.2;

    for i in 0..circle_s {
        let theta = (i as f64 / circle_s as f64) * 2.0 * PI + offset;
        let x = theta.cos();
        let y = ground + r;
        let z = theta.sin();

        scene.add(Sphere::new(
            DVec3::new(x, y, z),
            r,
            Material::specular(Texture::Solid(hsv_to_rgb(theta - offset)), 0.1),
        ));
    }

    let renderer = Renderer::new(scene, camera);
    renderer.render().save("circle.png")?;

    Ok(())
}
