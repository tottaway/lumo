use super::*;
use crate::{ Point, rand_utils };
use crate::tracer::ray::Ray;
use crate::tracer::object::{ Object, Disk };

/* testing framework to make sure materials don't give more energy than they receive,
 * "the white furnace"
 */

const NUM_RUNS: usize = 100;
const NUM_SAMPLES: usize = 10_000;
const MAX_RADIANCE: Float = 1.05;

fn white_texture() -> Texture {
    Texture::Solid(Color::WHITE)
}

fn disk(material: Material) -> Disk {
    *Disk::new(
        Point::ZERO,
        Normal::Z,
        1.0,
        material,
    )
}

#[test]
fn white_furnace_diffuse() {
    let m = Material::diffuse(white_texture());
    test_material(m);
}

#[test]
fn white_furnace_conductor75() {
    let m = Material::rough_mirror(white_texture(), 0.75, true);
    test_material(m);
}

#[test]
fn white_furnace_conductor50() {
    let m = Material::rough_mirror(white_texture(), 0.50, true);
    test_material(m);
}

#[test]
fn white_furnace_conductor25() {
    let m = Material::rough_mirror(white_texture(), 0.25, true);
    test_material(m);
}

#[test]
fn white_furnace_conductor10() {
    let m = Material::rough_mirror(white_texture(), 0.10, true);
    test_material(m);
}

#[test]
fn white_furnace_conductor0() {
    let m = Material::rough_mirror(white_texture(), 0.00, true);
    test_material(m);
}

#[test]
fn white_furnace_dielectric75_eta15() {
    let m = Material::transparent(white_texture(), 0.75, 1.5);
    test_material(m);
}

#[test]
fn white_furnace_dielectric50_eta15() {
    let m = Material::transparent(white_texture(), 0.50, 1.5);
    test_material(m);
}

#[test]
fn white_furnace_dielectric25_eta15() {
    let m = Material::transparent(white_texture(), 0.25, 1.5);
    test_material(m);
}

#[test]
fn white_furnace_dielectric10_eta15() {
    let m = Material::transparent(white_texture(), 0.10, 1.5);
    test_material(m);
}

#[test]
fn white_furnace_dielectric0_eta15() {
    let m = Material::transparent(white_texture(), 0.00, 1.5);
    test_material(m);
}

#[test]
fn white_furnace_dielectric75_eta25() {
    let m = Material::transparent(white_texture(), 0.75, 2.5);
    test_material(m);
}

#[test]
fn white_furnace_dielectric50_eta25() {
    let m = Material::transparent(white_texture(), 0.50, 2.5);
    test_material(m);
}

#[test]
fn white_furnace_dielectric25_eta25() {
    let m = Material::transparent(white_texture(), 0.25, 2.5);
    test_material(m);
}

#[test]
fn white_furnace_dielectric10_eta25() {
    let m = Material::transparent(white_texture(), 0.10, 2.5);
    test_material(m);
}

#[test]
fn white_furnace_dielectric0_eta25() {
    let m = Material::transparent(white_texture(), 0.00, 2.5);
    test_material(m);
}

fn test_material(m: Material) {
    let d = disk(m);
    for _ in 0..NUM_RUNS {
        assert!(furnace_pass(&d));
    }
}

fn furnace_pass(d: &Disk) -> bool {
    let origin = Point::Z;
    // random point on the disk
    let p = rand_utils::square_to_disk(rand_utils::unit_square());
    let r = Ray::new(origin, p.extend(0.0) - origin);
    let wo = r.dir;
    let radiance = match d.hit(&r, 0.0, 1e10) {
        None => Color::WHITE,
        Some(h) => furnace_sample(r, h),
    };
    let pass = radiance.rgb.max_element() < MAX_RADIANCE;

    if !pass {
        println!("L: {}, wo: {}", radiance.rgb, wo);
    }

    pass
}

fn furnace_sample(r: Ray, h: Hit) -> Color {
    let wo = r.dir;
    let m = h.material;

    let mut misses = 0;
    let mut radiance = Color::BLACK;

    for _ in 0..NUM_SAMPLES {
        match m.bsdf_sample(wo, &h, rand_utils::unit_square()) {
            None => misses += 1,
            Some(wi) => radiance += m.bsdf_f(wo, wi, Transport::Radiance, &h)
                * m.shading_cosine(wo, wi)
                / m.bsdf_pdf(wo, wi, &h, false),
        }
    }

    radiance / (NUM_SAMPLES - misses) as Float
}
