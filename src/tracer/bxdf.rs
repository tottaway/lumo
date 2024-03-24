use crate::{ Direction, Normal, Transport, Float, Vec2, rand_utils };
use crate::tracer::{ Color, microfacet::MfDistribution };

mod microfacet;
mod scatter;

#[derive(Clone, Copy)]
pub enum BxDF {
    Lambertian,
    /// Perfect transmission, to be merged with MfDielectric
    Transmission(Float),
    /// Perfect mirror, to be merged with MfReflection
    Reflection,
    /// Disney diffuse
    MfDiffuse(MfDistribution),
    /// Microfacet mirror
    MfReflection(MfDistribution),
    /// Microfacet glass
    MfDielectric(MfDistribution),
    None,
}

impl BxDF {
    pub fn is_specular(&self) -> bool {
        match self {
            Self::Reflection | Self::Transmission(_) | Self::MfDielectric(_) => true,
            Self::MfReflection(mfd) => mfd.is_specular(),
            _ => false,
        }
    }

    pub fn is_diffuse(&self) -> bool { !self.is_specular() }

    #[allow(clippy::match_like_matches_macro)]
    pub fn is_transmission(&self) -> bool {
        match self {
            Self::Transmission(_) | Self::MfDielectric(_) => true,
            _ => false
        }
    }

    pub fn is_reflection(&self) -> bool { !self.is_transmission() }

    pub fn is_delta(&self) -> bool {
        match self {
            Self::Reflection | Self::Transmission(_) => true,
            Self::MfReflection(mfd) | Self::MfDielectric(mfd) => mfd.is_delta(),
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
            Self::MfDiffuse(mfd) => microfacet::diffuse_f(wo, wi, mfd, albedo),
            Self::MfReflection(mfd) => microfacet::reflection_f(wo, wi, mfd, albedo),
            Self::MfDielectric(mfd) => microfacet::transmission_f(wo, wi, mfd, albedo, mode),
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
            Self::MfDielectric(mfd) => microfacet::transmission_sample(wo, mfd, rand_sq),
            Self::None => None,
        }
    }

    pub fn pdf(&self, wo: Direction, wi: Direction, swap_dir: bool) -> Float {
        match self {
            Self::Reflection | Self::Transmission(_) => 1.0,
            Self::Lambertian => scatter::lambertian_pdf(wi),
            Self::MfDiffuse(_) => scatter::lambertian_pdf(wi),
            Self::MfReflection(mfd) => microfacet::reflection_pdf(wo, wi, mfd),
            Self::MfDielectric(mfd) => microfacet::transmission_pdf(wo, wi, mfd, swap_dir),
            Self::None => 0.0,
        }
    }
}
