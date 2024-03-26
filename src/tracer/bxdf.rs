use crate::{ Direction, Normal, Transport, Float, Vec2, rand_utils };
use crate::tracer::{ Color, microfacet::MfDistribution };

mod microfacet;
mod scatter;

#[cfg(test)]
mod bxdf_tests;

#[cfg(test)]
mod bxdf_chi2_tests;

#[derive(Clone, Copy)]
pub enum BxDF {
    Lambertian,
    /// Lambertian diffuse
    MfDiffuse(MfDistribution),
    /// Microfacet mirror
    MfConductor(MfDistribution),
    /// Microfacet glass
    MfDielectric(MfDistribution),
    None,
}

impl BxDF {
    pub fn is_specular(&self) -> bool {
        match self {
            Self::MfDielectric(_) => true,
            Self::MfConductor(mfd) => mfd.is_specular(),
            _ => false,
        }
    }

    pub fn is_diffuse(&self) -> bool { !self.is_specular() }

    #[allow(clippy::match_like_matches_macro)]
    pub fn is_transmission(&self) -> bool {
        match self {
            Self::MfDielectric(_) => true,
            _ => false
        }
    }

    pub fn is_reflection(&self) -> bool { !self.is_transmission() }

    pub fn is_delta(&self) -> bool {
        match self {
            Self::MfConductor(mfd) | Self::MfDielectric(mfd) => mfd.is_delta(),
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
            Self::MfDiffuse(mfd) => microfacet::diffuse_f(wo, wi, mfd, albedo),
            Self::MfConductor(mfd) => microfacet::conductor_f(wo, wi, mfd, albedo),
            Self::MfDielectric(mfd) => microfacet::dielectric_f(wo, wi, mfd, albedo, mode),
            Self::None => Color::BLACK,
        }
    }

    pub fn sample(&self, wo: Direction, rand_sq: Vec2) -> Option<Direction> {
        match self {
            Self::Lambertian => Some( rand_utils::square_to_cos_hemisphere(rand_sq) ),
            Self::MfDiffuse(_) => Some( rand_utils::square_to_cos_hemisphere(rand_sq) ),
            Self::MfConductor(mfd) => microfacet::conductor_sample(wo, mfd, rand_sq),
            Self::MfDielectric(mfd) => microfacet::dielectric_sample(wo, mfd, rand_sq),
            Self::None => None,
        }
    }

    pub fn pdf(&self, wo: Direction, wi: Direction, swap_dir: bool) -> Float {
        match self {
            Self::Lambertian => scatter::lambertian_pdf(wi),
            Self::MfDiffuse(_) => scatter::lambertian_pdf(wi),
            Self::MfConductor(mfd) => microfacet::conductor_pdf(wo, wi, mfd),
            Self::MfDielectric(mfd) => microfacet::dielectric_pdf(wo, wi, mfd, swap_dir),
            Self::None => 0.0,
        }
    }
}
