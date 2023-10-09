use crate::{ Normal, Direction, Transport, Float, Vec3 };
use crate::tracer::{
    Color, hit::Hit, ray::Ray,
    microfacet::MfDistribution,
    texture::Texture, bsdf::BSDF, bxdf::BxDF,
};

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
    /// Metallic microfacet material
    pub fn metallic(texture: Texture, roughness: Float) -> Self {
        let mfd = MfDistribution::new(roughness, 1.5, 1.0);
        let bsdf = BSDF::new()
            .add(BxDF::MfDiffuse(mfd))
            .add(BxDF::MfReflection(mfd));
        Self::Microfacet(bsdf, texture, mfd)
    }

    /// Specular microfacet material
    pub fn specular(texture: Texture, roughness: Float) -> Self {
        let mfd = MfDistribution::new(roughness, 1.5, 0.0);
        let bsdf = BSDF::new()
            .add(BxDF::MfDiffuse(mfd))
            .add(BxDF::MfReflection(mfd));
        Self::Microfacet(bsdf, texture, mfd)
    }

    /// Diffuse material
    pub fn diffuse(texture: Texture) -> Self {
        let mfd = MfDistribution::new(1.0, 1.5, 0.0);
        let bsdf = BSDF::new()
            .add(BxDF::MfDiffuse(mfd))
            .add(BxDF::MfReflection(mfd));
        Self::Microfacet(bsdf, texture, mfd)
    }

    /// Transparent material
    pub fn transparent(texture: Texture, roughness: Float, eta: Float) -> Self {
        let mfd = MfDistribution::new(roughness, eta, 0.0);
        let bsdf = BSDF::new()
            .add(BxDF::MfTransmission(mfd));
        Self::Microfacet(bsdf, texture, mfd)
    }

    /// Perfect reflection
    pub fn mirror() -> Self {
        let bsdf = BSDF::new()
            .add(BxDF::Reflection);
        Self::Standard(bsdf, Texture::default())
    }

    /// Perfect refraction
    pub fn glass(eta: Float) -> Self {
        assert!(eta >= 1.0);
        let bsdf = BSDF::new()
            .add(BxDF::Transmission(eta));
        Self::Standard(bsdf, Texture::default())
    }

    /// Is the material specular? I.e. reflects light
    pub fn is_specular(&self) -> bool {
        false
        /*
        match self {
            Self::Mirror | Self::Glass(..) => true,
            Self::Microfacet(_, mfd) => mfd.is_specular(),
            _ => false,
    }
         */
    }

    /// Does the material scattering follow delta distribution?
    /// Dumb hack to make delta things not have shadows in path trace.
    pub fn is_delta(&self) -> bool {
        false
        /*
        match self {
            Self::Lambertian(_) => false,
            Self::Microfacet(_, mfd) => mfd.is_delta(),
            _ => true,
    }
         */
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

    // instead of this we want get BSDF, that updates normals albedos and such.
    // and gives REFERENCE to the updated BSDF (RACE CONDITIONS? yes it is a problem)
    // and how should we handle mediums then....
    // to avoid race conditions, maybe the BSDF itself should take normals and allbedo as args....
    /// What is the color at `h`?
    pub fn bsdf_f(
        &self,
        wo: Direction,
        wi: Direction,
        mode: Transport,
        h: &Hit
    ) -> Color {
        let ns = h.ns;
        let ng = h.ng;
        match self {
            Self::Standard(bsdf, texture) | Self::Microfacet(bsdf, texture, _) => {
                bsdf.update_normals(ns, ng);
                let albedo = texture.albedo_at(h);
                bsdf.f(wo, wi, albedo, mode)
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

    /// Computes the shading cosine coefficient per material
    pub fn shading_cosine(&self, wi: Direction, ns: Normal) -> Float {
        match self {
            //Self::Microfacet(..) | Self::Lambertian(_) => ns.dot(wi).abs(),
            _ => ns.dot(wi).abs()
        }
    }
}
