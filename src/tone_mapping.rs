use crate::DVec3;
use crate::rgb_to_luminance;

/// Enum for different tone mappers
pub enum ToneMap {
    /// Applies no tone mapping
    NoMap,
    /// Clamps values to \[0,1\]
    Clamp,
    /// Reinhard tone mapping to luminance `M(l) = l / (1 + l)`
    Reinhard,
    /// Hable tone map fron Uncharted 2
    HableFilmic,
    /// ACES (Academy Color Encoding System) approximation by Krzysztof Narkowicz
    ACES,
}

impl ToneMap {
    /// Tone maps the samples in `samples`
    pub fn map(&self, rgb: DVec3) -> DVec3 {
        match self {
            Self::NoMap => rgb,
            Self::Clamp => rgb.clamp(DVec3::ZERO, DVec3::ONE),
            Self::Reinhard => {
                let l_in = rgb_to_luminance(rgb);
                rgb / (1.0 + l_in)
            }
            Self::HableFilmic => {
                let exposure = 2.0;
                let curr = Self::hable_partial(rgb * exposure);
                let white = DVec3::splat(11.2);
                let white_scale = DVec3::ONE / Self::hable_partial(white);
                curr * white_scale
            }
            Self::ACES => {
                let rgb = rgb * 0.6;
                let a = 2.51;
                let b = 0.03;
                let c = 2.43;
                let d = 0.59;
                let e = 0.14;
                ((rgb * (a * rgb + b)) / (rgb * (c * rgb + d) + e))
                    .clamp(DVec3::ZERO, DVec3::ONE)
            }
        }
    }

    fn hable_partial(rgb: DVec3) -> DVec3 {
        let a = 0.15;
        let b = 0.50;
        let c = 0.10;
        let d = 0.20;
        let e = 0.02;
        let f = 0.30;

        (rgb * (a * rgb + c * b) + d * e)
            / (rgb * (a * rgb + b) + d * f)
            - e / f
    }
}