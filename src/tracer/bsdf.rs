use crate::{ Direction, Transport, Float, Vec2 };
use crate::tracer::{ Color, bxdf::BxDF, onb::Onb, hit::Hit };
use rand::prelude::SliceRandom;

#[allow(non_upper_case_globals)]
const MAX_BxDF: usize = 4;

#[allow(non_snake_case)]
pub struct BSDF {
    BxDFs: Vec<BxDF>,
}

impl BSDF {
    /// Construct new empty BSDF
    pub fn new() -> Self {
        Self {
            BxDFs: vec!(),
        }
    }

    /// Add `bxdf` to `self` (given we have less than `MAX_BxDF` BxDFs)
    pub fn add(mut self, bxdf: BxDF) -> Self {
        if self.BxDFs.len() < MAX_BxDF {
            self.BxDFs.push(bxdf);
        }
        self
    }

    pub fn is_specular(&self) -> bool {
        self.BxDFs.iter()
            .any(|bxdf| bxdf.is_specular())
    }

    /// Evaluate the BSDF
    pub fn f(
        &self,
        wo: Direction,
        wi: Direction,
        h: &Hit,
        albedo: Color,
        mode: Transport
    ) -> Color {
        let ns = h.ns;
        let uvw = Onb::new(ns);

        let wo_local = uvw.to_local(wo);
        let wi_local = uvw.to_local(wi);

        self.BxDFs.iter()
            .fold(Color::BLACK, |color, bxdf| {
                color + bxdf.f(wo_local, wi_local, albedo, mode)
            })
    }

    /// Sample direction from a random BxDF
    pub fn sample(
        &self,
        wo: Direction,
        h: &Hit,
        rand_sq: Vec2
    ) -> Option<Direction> {
        let ns = h.ns;
        let uvw = Onb::new(ns);

        let wo_local = uvw.to_local(wo);

        self.BxDFs.choose(&mut rand::thread_rng())
            .and_then(|bxdf| bxdf.sample(wo_local, rand_sq))
            .map(|wi| uvw.to_world(wi))
    }

    /// PDF for the BSDF
    pub fn pdf(
        &self,
        wo: Direction,
        wi: Direction,
        h: &Hit,
        swap_dir: bool
    ) -> Float {
        let ns = h.ns;
        let uvw = Onb::new(ns);

        let wo_local = uvw.to_local(wo);
        let wi_local = uvw.to_local(wi);

        self.BxDFs.iter()
            .fold(0.0, |prob, bxdf| {
                prob + bxdf.pdf(wo_local, wi_local, swap_dir)
            }) / self.BxDFs.len() as Float
    }
}
