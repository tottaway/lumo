use lumo::tracer::*;
use lumo::*;

fn main() -> Result<(), std::io::Error> {
    let camera = Camera::default(1024, 768);
    let def_color = Color::new(242, 242, 242);
    let mut scene = Scene::empty_box(
        def_color,
        // left
        Material::diffuse(Texture::Solid(Color::new(255, 0, 255))),
        // right
        Material::diffuse(Texture::Solid(Color::new(0, 255, 255))),
    );

    scene.add(Sphere::new(
        Vec3::new(-0.45, -0.5, -1.5),
        0.25,
        Material::mirror(),
    ));

    scene.add(Sphere::new(
        Vec3::new(0.45, -0.5, -1.3),
        0.25,
        Material::glass(2.5),
    ));

    let renderer = Renderer::new(scene, camera);
    renderer.render().save("box.png")?;
    Ok(())
}
