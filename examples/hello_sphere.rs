use rust_tracer::*;
use rust_tracer::tracer::*;
use glam::DVec3;

fn main() -> Result<(), png::EncodingError> {
    let camera = Camera::default();
    let mut scene = Scene::default();

    scene.add(
        Plane::new(
            DVec3::NEG_Y,
            DVec3::Y,
            Material::diffuse(Texture::Solid(srgb_to_linear(190, 200, 210)))
        )
    );

    scene.add(
        Sphere::new(
            8.0 * DVec3::Y + 1.5 * DVec3::NEG_Z,
            4.0,
            Material::Light(Texture::Solid(srgb_to_linear(255, 255, 255)))
        )
    );

    scene.add(
        Sphere::new(
            DVec3::ZERO,
            1.0,
            Material::diffuse(Texture::Solid(srgb_to_linear(0, 0, 255)))
        )
            .scale(0.3, 0.3, 0.3)
            .translate(0.0, -0.7, -1.5)
    );

    let mut renderer = Renderer::new(scene, camera);
    renderer.set_samples(36);
    renderer.render()
        .save("hello.png")
}
