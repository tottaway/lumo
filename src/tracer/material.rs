use crate::{ Normal, Direction, Transport, Float, Vec3, Vec2 };
use crate::tracer::{
    Color, hit::Hit, microfacet::MfDistribution,
    texture::Texture, bsdf::BSDF, bxdf::BxDF, onb::Onb
};

#[cfg(test)]
mod white_furnace_tests;

/// Describes which material an object is made out of
pub enum Material {
    /// Material with microfacet BxDF(s)
    Microfacet(BSDF, Texture, MfDistribution),
    /// Material without microfacet BxDF(s)
    Standard(BSDF, Texture),
    /// Emits light
    Light(Texture),
    /// Volumetric material for mediums. `scatter_param`, `sigma_t`, `sigma_s`
    Volumetric(Float, Vec3, Color),
    /// Not specified. Used with objects that are built on top of other objects.
    Blank,
}

impl Material {
    /// General microfacet constructor
    pub fn microfacet(
        texture: Texture,
        roughness: Float,
        eta: Float,
        k: Float,
        is_transparent: bool,
        fresnel_enabled: bool
    ) -> Self {
        let mfd = MfDistribution::new(roughness, eta, k, fresnel_enabled);
        // dirty dirty...
        let bsdf = if is_transparent {
            BSDF::new(BxDF::MfDielectric(mfd))
        } else if fresnel_enabled {
            BSDF::new(BxDF::MfConductor(mfd))
        } else {
            BSDF::new(BxDF::MfDiffuse(mfd))
        };
        Self::Microfacet(bsdf, texture, mfd)
    }

    /// Microfacet mirror with assignable roughness
    pub fn metal(texture: Texture, roughness: Float, eta: Float, k: Float) -> Self {
        let is_transparent = false;
        let fresnel_enabled = true;

        Self::microfacet(
            texture,
            roughness,
            eta,
            k,
            is_transparent,
            fresnel_enabled,
        )
    }

    /// Diffuse material
    pub fn diffuse(texture: Texture) -> Self {
        let roughness = 1.0;
        let eta = 1.5;
        let k = 0.0;
        let is_transparent = false;
        let fresnel_enabled = false;

        Self::microfacet(
            texture,
            roughness,
            eta,
            k,
            is_transparent,
            fresnel_enabled,
        )
    }

    /// Transparent material
    pub fn transparent(texture: Texture, roughness: Float, eta: Float) -> Self {
        let k = 0.0;
        let is_transparent = true;
        let fresnel_enabled = true;

        Self::microfacet(
            texture,
            roughness,
            eta,
            k,
            is_transparent,
            fresnel_enabled,
        )
    }

    /// Perfect reflection
    pub fn mirror() -> Self {
        let texture = Texture::Solid(Color::WHITE);
        let roughness = 0.0;
        let eta = 1e5;
        let k = 0.0;
        let is_transparent = false;
        let fresnel_enabled = true;

        Self::microfacet(
            texture,
            roughness,
            eta,
            k,
            is_transparent,
            fresnel_enabled,
        )
    }

    /// Perfect refraction
    pub fn glass(eta: Float) -> Self {
        let texture = Texture::Solid(Color::WHITE);
        let roughness = 0.0;
        let k = 0.0;
        let is_transparent = true;
        let fresnel_enabled = true;

        Self::microfacet(
            texture,
            roughness,
            eta,
            k,
            is_transparent,
            fresnel_enabled,
        )
    }

    /// Is the material specular? I.e. reflects light
    pub fn is_specular(&self) -> bool {
        match self {
            Self::Volumetric(..) => true,
            Self::Microfacet(bsdf, _, _) | Self::Standard(bsdf, _) => bsdf.is_specular(),
            _ => false,
        }
    }

    /// Does the material scattering follow delta distribution?
    /// Dumb hack to make delta things not have shadows in path trace.
    pub fn is_delta(&self) -> bool {
        match self {
            Self::Microfacet(bsdf, _, _) | Self::Standard(bsdf, _) => bsdf.is_delta(),
            _ => false,
        }
    }


    /// How much light emitted at `h`?
    pub fn emit(&self, h: &Hit) -> Color {
        match self {
            Self::Light(t) => if h.backface {
                Color::BLACK
            } else {
                t.albedo_at(h)
            },
            _ => Color::BLACK
        }
    }

    /// BSDF evaluated at `h` for incoming `wo` and outgoing `wi` while
    /// transporting `mode`
    pub fn bsdf_f(
        &self,
        wo: Direction,
        wi: Direction,
        mode: Transport,
        h: &Hit
    ) -> Color {
        match self {
            Self::Standard(bsdf, texture) | Self::Microfacet(bsdf, texture, _) => {
                let albedo = texture.albedo_at(h);
                bsdf.f(wo, wi, h, albedo, mode)
            }
            // volumetric BSDF handled in integrator to cancel out PDF
            Self::Volumetric(_, sigma_t, sigma_s) => {
                let transmittance = (-*sigma_t * h.t).exp();
                // cancel out the transmittance pdf taken from scene transmitance
                let pdf = (transmittance * *sigma_t).dot(Vec3::ONE)
                    / transmittance.dot(Vec3::ONE);

                if pdf == 0.0 { Color::WHITE } else { *sigma_s / pdf }
            }
            _ => Color::BLACK,
        }
    }

    /// Samples leaving direction from `h` from incoming direction `wo`
    pub fn bsdf_sample(
        &self,
        wo: Direction,
        h: &Hit,
        rand_sq: Vec2
    ) -> Option<Direction> {
        match self {
            Self::Standard(bsdf, _) | Self::Microfacet(bsdf, _, _) => {
                bsdf.sample(wo, h, rand_sq)
            }
            /* Henyey-Greenstein (1941) */
            Self::Volumetric(g, _, _) => {
                let cos_theta = if g.abs() < 1e-3 {
                    1.0 - 2.0 * rand_sq.x
                } else {
                    let g2 = g * g;
                    let fract = (1.0 - g2)
                        / (1.0 - g + 2.0 * g * rand_sq.x);
                    (1.0 + g2 - fract * fract) / (2.0 * g)
                };
                let sin_theta = (1.0 - cos_theta * cos_theta).max(0.0).sqrt();

                let phi = 2.0 * crate::PI * rand_sq.y;

                let v = -wo;
                let uvw = Onb::new(v);
                let wi = uvw.to_world(Direction::new(
                    sin_theta * phi.cos(),
                    sin_theta * phi.sin(),
                    cos_theta
                ));

                Some( wi )
            }
            _ => None,
        }
    }

    /// PDF for direction `wi` at `h` with incoming direction `wo`
    pub fn bsdf_pdf(
        &self,
        wo: Direction,
        wi: Direction,
        h: &Hit,
        swap_dir: bool
    ) -> Float {
        match self {
            Self::Standard(bsdf, _) | Self::Microfacet(bsdf, _, _) => {
                bsdf.pdf(wo, wi, h, swap_dir)
            }
            Self::Volumetric(g, _, _) => {
                let v = -wo;
                let cos_theta = v.dot(wi);

                let g2 = g * g;
                let denom = 1.0 + g2 + 2.0 * g * cos_theta;

                (1.0 - g2) / (4.0 * crate::PI * denom * denom.max(0.0).sqrt())
            }
            _ => 0.0,
        }
    }

    /// Computes the shading cosine coefficient per material
    pub fn shading_cosine(&self, wi: Direction, ns: Normal) -> Float {
        match self {
            Self::Microfacet(..) => ns.dot(wi).abs(),
            _ => 1.0,
        }
    }
}
