use crate::{ Direction, Point, Normal, Transport, Float, Vec2 };
use crate::tracer::{ Color, ray::Ray, bxdf::BxDF };

pub struct BSDF {
    BxDFs: Vec<BxDF>,
}

impl BSDF {
    pub fn new() -> Self {
        Self { BxDFs: vec![] }
    }

    pub fn add(&mut self, bxdf: BxDF) {
        self.BxDFs.push(bxdf)
    }

    pub fn eval(&self, wo: Direction, wi: Direction, mode: Transport) -> Color {
        // local transform
        self.BxDFs.iter()
            .fold(Color::BLACK, |color, bxdf| color + bxdf.eval(wo, wi, mode))
    }

    pub fn sample(&self, wo: Direction, rand_sq: Vec2) -> Option<Direction> {
        self.BxDFs.choose(rand::thread_rng())
            .and_then(|bxdf| bxdf.sample(wo, rand_sq))
    }

    pub fn pdf(&self, wo: Direction, wi: Direction, swap_dir: bool) -> Float {
        self.BxDFs.iter()
            .fold(0.0, |prob, bxdf| prob + bxdf.pdf(wo, wi, swap_dir))
            / self.BxDFs.len()
    }
}
