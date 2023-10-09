use crate::{ Direction, Normal, Transport, Float, Vec2, rand_utils };
use crate::tracer::{ Color, ray::Ray, microfacet::MfDistribution };

mod microfacet;
mod scatter;

pub enum BxDF {
    Lambertian,
    Transmission(Float),
    Reflection,
    MfReflection(MfDistribution),
    MfTransmission(MfDistribution),
}

impl BxDF {
    pub fn f(&self, wo: Direction, wi: Direction, albedo: Color, mode: Transport) -> Color {
        match self {
            Self::Lambertian => albedo / crate::PI,
            Self::Reflection => Color::WHITE,
            Self::Transmission(eta) => scatter::transmission_f(wo, *eta, mode),
            Self::MfReflection(mfd) => microfacet::reflection_f(wo, wi, mfd, albedo),
            Self::MfTransmission(mfd) => microfacet::transmission_f(wo, wi, mfd, albedo, mode),
        }
    }

    pub fn sample(&self, wo: Direction, rand_sq: Vec2) -> Option<Direction> {
        match self {
            Self::Reflection => Some( Direction::new(wo.x, wo.y, -wo.z) ),
            Self::Lambertian => Some( rand_utils::square_to_cos_hemisphere(rand_sq) ),
            Self::Transmission(eta) => scatter::transmission_sample(wo, *eta, rand_sq),
        }
    }

    pub fn pdf(&self, wo: Direction, wi: Direction, swap_dir: bool) -> Float {
        match self {
            Self::Reflection | Self::Transmission(_) => 1.0,
            Self::Lambertian => scatter::lambertian_pdf(wi),
        }
    }
}
