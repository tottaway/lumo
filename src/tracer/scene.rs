use crate::{DVec3, DMat3, DAffine3, DVec2};
use std::f64::consts::PI;
#[allow(unused_imports)]
use crate::samplers::{UniformSampler, JitteredSampler};
use crate::perlin::Perlin;
use crate::consts::{EPSILON, SHADOW_RAYS};
use crate::tracer::hit::Hit;
use crate::tracer::ray::Ray;
use crate::tracer::texture::Texture;
use crate::tracer::material::Material;
use crate::tracer::object::sphere::Sphere;
use crate::tracer::object::{Object, Plane, Rectangle, Cuboid};

#[cfg(test)]
mod scene_tests;

pub struct Scene {
    pub ambient: DVec3,
    /* vec of indices to objects that are lights */
    lights: Vec<usize>,
    /* TODO */
    pub objects: Vec<Box<dyn Object>>,
}

/* temporary constant */
const LIGHT_R: f64 = 0.1;

type PixelSampler = JitteredSampler;

impl Scene {
    pub fn new(amb: DVec3, objs: Vec<Box<dyn Object>>) -> Self {
        let lights = (0..objs.len()).map(|i: usize| match objs[i].material() {
            Material::Light(_) => i,
            _ => objs.len(),
        }).filter(|i: &usize| *i != objs.len()).collect();

        Self {
            ambient: amb,
            lights: lights,
            objects: objs,
        }
    }

    /* might want to print x of this, y of that, ... */
    pub fn size(&self) -> usize {
        self.objects.iter().fold(0, |acc, obj| acc + obj.size())
    }

    pub fn hit(&self, r: &Ray) -> Option<Hit> {
        self.objects.iter().map(|obj| obj.hit(r))
            .fold(None, |closest, hit| {
                if closest.is_none() || (hit.is_some() && hit < closest) {
                    hit
                } else {
                    closest
                }
            })
    }

    pub fn rays_to_light(&self, h: &Hit) -> Vec<Ray> {
        self.lights.iter().flat_map(|light_idx: &usize| {
            PixelSampler::new(SHADOW_RAYS).map(|p: DVec2| {
                self.objects[*light_idx].sample_ray(h, p)
            }).filter(|r: &Ray| self.hit_light(r, &self.objects[*light_idx]))
                .collect::<Vec<Ray>>()
        }).collect()
    }

    fn hit_light(&self, r: &Ray, light: &Box<dyn Object>) -> bool {
        let lt = light.hit(r).map_or(f64::INFINITY, |h| h.t);
        let no_block_light = |obj: &&Box<dyn Object>| -> bool {
            /* only keep hits that are closer than light and not translucent */
            obj.hit(r).filter(|h: &Hit| {
                h.t + EPSILON < lt && !h.object.is_translucent()
                /* to check if object is behind light the rays are same,
                 * so we can just compare t's. for our light the t is
                 * always 1.0 by design of r */
            }).is_none()
        };

        self.objects.iter().take_while(no_block_light)
            .count() == self.objects.len()
        }

        pub fn box_scene() -> Self {
            /* y ground */
            let yg = -0.8;
            let col = DVec3::new(255.0, 253.0, 208.0) / 255.9;
            let light_xy = 0.2;
            let light_z = yg - light_xy;
            Self::new(
                DVec3::splat(0.0),
                vec![
                    /* if light and roof have same y, then light should come
                     * before roof in this vector and everything should be ok */
                    Rectangle::new(
                        DMat3::from_cols(
                            DVec3::new(-light_xy, -yg, light_z - light_xy),
                            DVec3::new(-light_xy, -yg, light_z + light_xy),
                            DVec3::new(light_xy, -yg, light_z + light_xy),
                        ),
                        Material::Light(Texture::Solid(DVec3::ONE)),
                    ),
                    Sphere::new(
                        DVec3::new(-0.4, -0.6, -1.2),
                        0.2,
                        Material::Mirror,
                    ),
                    Cuboid::new(
                        DAffine3::from_translation(
                            DVec3::new(0.2, yg, 1.7*yg))
                            * DAffine3::from_scale(DVec3::new(0.2, 0.4, 0.2))
                            * DAffine3::from_rotation_y(PI / 10.0),
                        Material::Phong(
                            Texture::Solid(DVec3::new(0.0, 0.9, 0.0))
                        )
                    ),
                    // roof
                    Plane::new(
                        DVec3::new(0.0, -yg, 0.0),
                        DVec3::new(0.0, -1.0, 0.0),
                        Material::Phong(Texture::Solid(col)),
                    ),
                    /* floor */
                    Rectangle::new(
                        DMat3::from_cols(
                            DVec3::new(-yg, yg, 0.0),
                            DVec3::new(yg, yg, 0.0),
                            DVec3::new(yg, yg, 2.0*yg),
                        ),
                        Material::Phong(Texture::Solid(col)),
                    ),
                    // front wall
                    Plane::new(
                        DVec3::new(0.0, 0.0, 2.0*yg),
                        DVec3::new(0.0, 0.0, 1.0),
                        Material::Phong(Texture::Solid(col)),
                    ),
                    // left wall
                    Plane::new(
                        DVec3::new(yg, 0.0, 0.0),
                        DVec3::new(1.0, 0.0, 0.0),
                        Material::Phong(Texture::Solid(DVec3::new(0.0, 1.0, 1.0))),
                    ),
                    // right wall
                    Plane::new(
                        DVec3::new(-yg, 0.0, 0.0),
                        DVec3::new(-1.0, 0.0, 0.0),
                        Material::Phong(Texture::Solid(DVec3::new(1.0, 0.0, 1.0))),
                    ),
                    // background
                    Plane::new(
                        DVec3::new(0.0, 0.0, 0.1),
                        DVec3::new(0.0, 0.0, -1.0),
                        Material::Blank,
                    ),
                ],
            )
        }

        pub fn default() -> Self {
            Self::new(
                DVec3::splat(0.0),
                vec![
                    Sphere::new(
                        DVec3::new(-0.3, 0.2, -0.1),
                        LIGHT_R,
                        Material::Light(Texture::Solid(DVec3::ONE)),
                    ),
                    // floor
                    Plane::new(
                        DVec3::new(0.0, -0.5, 0.0),
                        DVec3::new(0.0, 1.0, 0.0),
                        Material::Phong(Texture::Checkerboard(
                            Box::new(Texture::Checkerboard(
                                Box::new(Texture::Solid(DVec3::ZERO)),
                                Box::new(Texture::Solid(DVec3::ONE)),
                                4.0,
                            )),
                            /* share same perlin between all textures?
                             * could make cool checkers that way */
                            Box::new(Texture::Marble(
                                Perlin::new(DVec3::splat(192.0) / 255.9)
                            )),
                            1.0,
                        )),
                    ),
                    // right
                    Plane::new(
                        DVec3::new(3.0, 0.0, -3.0),
                        DVec3::new(-1.0, 0.0, 1.0),
                        Material::Phong(Texture::Solid(DVec3::new(0.0, 0.0, 1.0))),
                    ),
                    Rectangle::new(
                        DMat3::from_cols(
                            DVec3::new(1.2, 0.2, -0.8),
                            DVec3::new(0.8, 0.6, -0.4),
                            DVec3::new(0.4, 0.6, -0.8),
                        ),
                        Material::Mirror,
                    ),
                    Cuboid::new(
                        DAffine3::from_translation(
                            DVec3::new(-0.4, -0.5, -0.6))
                            * DAffine3::from_scale(DVec3::splat(0.15))
                            * DAffine3::from_rotation_y(PI / 4.0),
                        Material::Phong(Texture::Checkerboard(
                            Box::new(Texture::Solid(DVec3::new(1.0, 0.0, 1.0))),
                            Box::new(Texture::Solid(
                                DVec3::new(50.0, 205.0, 50.0) / 255.9
                            )),
                            9.0,
                        )),
                    ),
                    // left
                    Plane::new(
                        DVec3::new(-3.0, 0.0, -3.0),
                        DVec3::new(1.0, 0.0, 1.0),
                        Material::Phong(Texture::Solid(DVec3::new(1.0, 0.0, 0.0))),
                    ),
                    // behind
                    Plane::new(
                        DVec3::new(0.0, 0.0, 1.0),
                        DVec3::new(0.0, 0.0, -1.0),
                        Material::Phong(Texture::Solid(DVec3::new(1.0, 0.0, 1.0))),
                    ),
                    Sphere::new(
                        DVec3::new(0.0, 0.0, -1.0),
                        0.5,
                        Material::Phong(Texture::Solid(
                            DVec3::new(136.0, 8.0, 8.0) / 255.9
                        )),
                    ),
                    Sphere::new(
                        DVec3::new(-0.9, 0.0, -1.0),
                        0.1,
                        Material::Mirror,
                    ),
                    Sphere::new(
                        DVec3::new(-0.4, -0.12, -0.5),
                        0.1,
                        Material::Glass,
                    ),
                    Sphere::new(
                        DVec3::new(0.4, -0.2, -0.5),
                        0.1,
                        Material::Phong(Texture::Marble(Perlin::new(
                            DVec3::new(255.0, 182.0, 193.0) / 255.9
                        ))),
                    ),
                ]
            )
        }
    }
