use crate::{ Direction, Transport, Float, Vec2 };
use crate::tracer::{ Color, bxdf::BxDF, onb::Onb, hit::Hit };

#[allow(non_snake_case)]
pub struct BSDF {
    BxDF: BxDF
}

impl BSDF {
    /// Construct new empty BSDF
    #[allow(non_snake_case)]
    pub fn new(BxDF: BxDF) -> Self {
        Self { BxDF }
    }

    pub fn is_specular(&self) -> bool {
        self.BxDF.is_specular()
    }

    pub fn is_delta(&self) -> bool {
        self.BxDF.is_delta()
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

        self.BxDF.f(wo_local, wi_local, albedo, mode)
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

        self.BxDF.sample(wo_local, rand_sq)
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

        self.BxDF.pdf(wo_local, wi_local, swap_dir)
    }
}
