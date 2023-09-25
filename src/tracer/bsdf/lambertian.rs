use super::*;

/// Lambertian diffuse material
pub struct Lambertian {
    /// Albedo at the point of impact
    albedo: Color,
}

impl Glass {
    pub fn new(ro: &Ray, albedo: Color, ) -> Self {
        let v = -ro.dir;
        let ng = ho.ng;
        let inside = ng.dot(v) < 0.0;
        let eta_ratio = if inside { eta } else { 1.0 / eta };
        let ns = if inside { -ho.ns } else { ho.ns };

        let wi = refract(eta_ratio, v, ns);

        Self { wi, eta, ng }
    }
}

impl BSDF for Glass {
    fn eval(&self, _wo: Direction, wi: Direction, mode: Transport) -> Color {
        match mode {
            Transport::Importance => Color::WHITE,
            Transport::Radiance => {
                let inside = wi.dot(self.ng) > 0.0;
                if inside {
                    Color::splat(1.0 / (eta * eta))
                } else {
                    Color::splat(eta * eta)
                }
            }
        }
    }

    fn sample_direction(&self, _rand_sq: Vec2) -> Option<Direction> {
        Some( self.wi )
    }

    fn pdf_for(&self, ri: &Ray, swap_dir: bool) -> Float {
        let wi = ri.dir;
        if wi.dot(self.wi) >= 1.0 - crate::EPSILON { 1.0 } else { 0.0 }
    }
}
