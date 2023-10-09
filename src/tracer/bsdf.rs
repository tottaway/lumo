use crate::{ Direction, Point, Normal, Transport, Float, Vec2 };
use crate::tracer::{ Color, ray::Ray, bxdf::BxDF, onb::Onb };
use rand::prelude::SliceRandom;

pub struct BSDF {
    BxDFs: Vec<BxDF>,
    uvw: Onb,
    ng: Normal,
}

impl BSDF {
    pub fn new() -> Self {
        Self {
            BxDFs: vec![],
            uvw: Onb::new(Normal::Z),
            ng: Normal::Z,
        }
    }

    pub fn add(&mut self, bxdf: BxDF) {
        self.BxDFs.push(bxdf)
    }

    pub fn update_normals(&mut self, ns: Normal, ng: Normal) {
        self.uvw = Onb::new(ns);
        self.ng = ng;
    }

    pub fn eval(&self, wo: Direction, wi: Direction, mode: Transport) -> Color {
        let wo_local = self.uvw.to_local(wo);
        let wi_local = self.uvw.to_local(wi);

        self.BxDFs.iter()
            .fold(Color::BLACK, |color, bxdf| {
                color + bxdf.f(wo_local, wi_local, Color::BLACK, mode)
            })
    }

    pub fn sample(&self, wo: Direction, rand_sq: Vec2) -> Option<Direction> {
        let wo_local = self.uvw.to_local(wo);

        self.BxDFs.choose(&mut rand::thread_rng())
            .and_then(|bxdf| bxdf.sample(wo_local, rand_sq))
    }

    pub fn pdf(&self, wo: Direction, wi: Direction, swap_dir: bool) -> Float {
        let wo_local = self.uvw.to_local(wo);
        let wi_local = self.uvw.to_local(wi);

        self.BxDFs.iter()
            .fold(0.0, |prob, bxdf| {
                prob + bxdf.pdf(wo_local, wi_local, swap_dir)
            }) / self.BxDFs.len() as Float
    }
}
