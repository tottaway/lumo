use crate::{ Direction, Normal, Transport, Float, Vec2, rand_utils };
use crate::tracer::{ Color, microfacet::MfDistribution };

mod microfacet;
mod scatter;

#[derive(Clone, Copy)]
pub enum BxDF {
    Lambertian,
    Transmission(Float),
    Reflection,
    MfDiffuse(MfDistribution),
    MfReflection(MfDistribution),
    MfTransmission(MfDistribution),
    None,
}

impl BxDF {
    pub fn is_specular(&self) -> bool {
        match self {
            Self::Reflection | Self::Transmission(_) | Self::MfTransmission(_) => true,
            Self::MfReflection(mfd) => mfd.is_specular(),
            _ => false,
        }
    }

    pub fn f(
        &self,
        wo: Direction,
        wi: Direction,
        albedo: Color,
        mode: Transport
    ) -> Color {
        match self {
            Self::Lambertian => albedo / crate::PI,
            Self::Reflection => Color::WHITE,
            Self::Transmission(eta) => scatter::transmission_f(wo, *eta, mode),
            Self::MfDiffuse(mfd) => microfacet::diffuse_f(wo, wi, albedo, mfd),
            Self::MfReflection(mfd) => microfacet::reflection_f(wo, wi, mfd, albedo),
            Self::MfTransmission(mfd) => microfacet::transmission_f(wo, wi, mfd, albedo, mode),
            Self::None => Color::BLACK,
        }
    }

    pub fn sample(&self, wo: Direction, rand_sq: Vec2) -> Option<Direction> {
        match self {
            Self::Reflection => Some( Direction::new(wo.x, wo.y, -wo.z) ),
            Self::Lambertian => Some( rand_utils::square_to_cos_hemisphere(rand_sq) ),
            Self::Transmission(eta) => scatter::transmission_sample(wo, *eta),
            Self::MfDiffuse(_) => Some( rand_utils::square_to_cos_hemisphere(rand_sq) ),
            Self::MfReflection(mfd) => microfacet::reflection_sample(wo, mfd, rand_sq),
            Self::MfTransmission(mfd) => microfacet::transmission_sample(wo, mfd, rand_sq),
            Self::None => None,
        }
    }

    pub fn pdf(&self, wo: Direction, wi: Direction, swap_dir: bool) -> Float {
        match self {
            Self::Reflection | Self::Transmission(_) => 1.0,
            Self::Lambertian => scatter::lambertian_pdf(wi),
            Self::MfDiffuse(_) => scatter::lambertian_pdf(wi),
            Self::MfReflection(mfd) => microfacet::reflection_pdf(wo, wi, mfd),
            Self::MfTransmission(mfd) => microfacet::transmission_pdf(wo, wi, mfd, swap_dir),
            Self::None => 0.0,
        }
    }
}
