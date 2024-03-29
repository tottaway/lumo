use crate::tracer::onb::Onb;
use crate::{ Normal, Direction, Float, Vec2, spherical_utils };
use num::complex::Complex;

/// Configurable parameters for a microsurface
#[derive(Copy, Clone)]
pub struct MicrofacetConfig {
    /// Roughness of the surface (α) [0,1]
    pub roughness: Float,
    /// Refraction index of the material >= 1.0
    pub eta: Float,
    /// Absoprtion coefficient
    pub k: Float,
    /// Is fresnel term enabled
    pub fresnel_enabled: bool,
}

impl MicrofacetConfig {
    pub fn new(
        roughness: Float,
        eta: Float,
        k: Float,
        fresnel_enabled: bool,
    ) -> Self {
        assert!((0.0..=1.0).contains(&roughness));
        //assert!((0.0..=1.0).contains(&metallicity));
        assert!(eta >= 1.0);

        Self {
            roughness: roughness.max(1e-5),
            eta,
            k,
            fresnel_enabled,
        }
    }
}

/// Defines a distribution of normals for a microfacet. `Float` parameter is the
/// roughness (α) of the surface.
#[derive(Copy, Clone)]
pub enum MfDistribution {
    /// Walter et al. 2007
    Ggx(MicrofacetConfig),
    /// Beckmann et al. 1987
    Beckmann(MicrofacetConfig),
}

impl MfDistribution {
    pub fn new(
        roughness: Float,
        eta: Float,
        metallicity: Float,
        fresnel_enabled: bool
    ) -> Self {
        Self::Ggx(MicrofacetConfig::new(roughness, eta, metallicity, fresnel_enabled))
    }

    /// might need tuning, send ratio that emittance is multiplied with?
    pub fn is_specular(&self) -> bool {
        self.get_config().roughness < 0.01
    }

    /// Does the material have delta scattering distribution?
    pub fn is_delta(&self) -> bool {
        self.get_config().roughness < 1e-3
    }

    /// Get refraction index from config
    pub fn eta(&self) -> Float {
        self.get_config().eta
    }

    /// Get absorption coefficient from config
    pub fn k(&self) -> Float {
        self.get_config().k
    }

    /// Get roughness from config
    pub fn roughness(&self) -> Float {
        self.get_config().roughness
    }

    pub fn fresnel_enabled(&self) -> bool {
        self.get_config().fresnel_enabled
    }

    /// Getter, better way to do this?
    fn get_config(&self) -> &MicrofacetConfig {
        match self {
            Self::Ggx(cfg) | Self::Beckmann(cfg) => cfg,
        }
    }

    /// Disney diffuse (Burley 2012) with renormalization to conserve energy
    /// as done in Frostbite (Lagarde et al. 2014)
    pub fn _disney_diffuse(
        &self,
        cos_theta_v: Float,
        cos_theta_wh: Float,
        cos_theta_wi: Float
    ) -> Float {
        let roughness2 = self.roughness().powi(2);
        let energy_bias = 0.5 * roughness2;
        let fd90 = energy_bias + 2.0 * cos_theta_wh.powi(2) * roughness2;

        let view_scatter = 1.0 + (fd90 - 1.0) * (1.0 - cos_theta_v).powi(5);
        let light_scatter = 1.0 + (fd90 - 1.0) * (1.0 - cos_theta_wi).powi(5);

        let energy_factor = 1.0 + roughness2 * (1.0 / 1.51 - 1.0);

        view_scatter * light_scatter * energy_factor
    }

    /// The microfacet distribution function.
    ///
    /// # Distributions
    /// * Beckmann - exp(-tan^2(θ) / α^2) / (π * α^2 * cos^4(θ))
    /// * GGX - α^2 / (π * (cos^4(θ) * (α^2 - 1.0) + 1.0)^2)
    ///
    /// # Arguments
    /// * `wh` - Microsurface normal in shading space
    pub fn d(&self, wh: Normal) -> Float {
        match self {
            Self::Ggx(cfg) => {
                let tan2_theta = spherical_utils::tan2_theta(wh);

                if tan2_theta.is_nan() {
                    0.0
                } else {
                    let cos4_theta = spherical_utils::cos2_theta(wh).powi(2);
                    if cos4_theta < crate::EPSILON {
                        return 0.0;
                    }
                    let cos_phi = spherical_utils::cos_phi(wh);
                    let sin_phi = spherical_utils::sin_phi(wh);

                    let alpha2 = (cfg.roughness * cos_phi).powi(2)
                        + (cfg.roughness * sin_phi).powi(2);

                    alpha2 / (crate::PI * cos4_theta * (alpha2 + tan2_theta).powi(2))
                }
            }
            Self::Beckmann(cfg) => {
                let tan2_theta = spherical_utils::tan2_theta(wh);

                if tan2_theta.is_nan() {
                    0.0
                } else {
                    let roughness2 = cfg.roughness * cfg.roughness;
                    let cos4_theta = spherical_utils::cos2_theta(wh).powi(2);

                    (-tan2_theta / roughness2).exp()
                        / (crate::PI * roughness2 * cos4_theta)
                }
            }
        }
    }

    /// Fresnel term with the full equations
    /// # Arguments
    /// * `v`      - Direction to viewer in shading space
    /// * `wh`     - Microsurface normal in shading space
    pub fn f(&self, v: Direction, wh: Normal) -> Float {
        if self.k() == 0.0 {
            // use simpler form of fresnel equations for dielectrics
            self.fr_transmission(v, wh)
        } else {
            self.fr_complex(v, wh)
        }
    }

    fn fr_complex(&self, v: Direction, wh: Normal) -> Float {
        // this is a complex number: n + ik
        let eta = Complex::new(self.eta(), self.k());
        let cos_o = v.dot(wh).clamp(0.0, 1.0);
        let sin2_o = 1.0 - cos_o * cos_o;

        let sin2_i: Complex<Float> = sin2_o / (eta * eta);
        let cos_i: Complex<Float> = (1.0 - sin2_i).sqrt();

        let r_par: Complex<Float> = (eta * cos_o - cos_i) / (eta * cos_o + cos_i);
        let r_per: Complex<Float> = (cos_o - eta * cos_i) / (cos_o + eta * cos_i);

        (r_par.norm() + r_per.norm()) / 2.0
    }

    fn fr_transmission(&self, v: Direction, wh: Normal) -> Float {
        let cos_o = v.dot(wh);
        let inside = cos_o < 0.0;
        let (eta_o, eta_i) = if inside {
            (self.eta(), 1.0)
        } else {
            (1.0, self.eta())
        };
        let cos_o = cos_o.abs();
        let sin_o = (1.0 - cos_o * cos_o).max(0.0).sqrt();
        let sin_i = eta_o / eta_i * sin_o;

        // total internal reflection
        if sin_i >= 1.0 {
            return 1.0;
        }

        let cos_i = (1.0 - sin_i * sin_i).max(0.0).sqrt();

        let r_par = ( (eta_i * cos_o) - (eta_o * cos_i) )
                  / ( (eta_i * cos_o) + (eta_o * cos_i) );
        let r_per = ( (eta_o * cos_o) - (eta_i * cos_i) )
                  / ( (eta_o * cos_o) + (eta_i * cos_i) );

        (r_par * r_par + r_per * r_per) / 2.0
    }

    fn chi_pass(&self, v: Direction, wh: Normal) -> bool {
        // signum to fix refraction
        let cos_theta_wh = spherical_utils::cos_theta(wh);
        let cos_theta_v = spherical_utils::cos_theta(v);
        let chi = cos_theta_wh.signum() * v.dot(wh) * cos_theta_v;
        chi > crate::EPSILON
    }


    /// Shadow-masking term. Used to make sure that only microfacets that are
    /// visible from `v` direction are considered. Uses the method described
    /// in Chapter 8.4.3 of PBR due to Heitz et al. 2013.
    ///
    /// # Arguments
    /// * `v`  - View direction in shading space
    /// * `wi` - Direction of ray away from the point of impact in shading space
    /// * `wh` - Microsurface normal in shading space
    pub fn g(&self, v: Direction, wi: Direction, wh: Normal) -> Float {
        if !self.chi_pass(v, wh) {
            0.0
        } else {
            1.0 / (1.0 + self.lambda(v) + self.lambda(wi))
        }
    }

    pub fn g1(&self, v: Direction, wh: Normal) -> Float {
        if !self.chi_pass(v, wh) {
            0.0
        } else {
            1.0 / (1.0 + self.lambda(v))
        }
    }

    /// Lambda function used in the definition of the shadow-masking term.
    /// Beckmann with polynomial approximation and GGX exactly. PBR Chapter 8.4.3
    ///
    /// # Arguments
    /// * `w` - Direction to consider in shading space
    fn lambda(&self, w: Direction) -> Float {
        match self {
            Self::Ggx(cfg) => {
                let tan2_theta = spherical_utils::tan2_theta(w);

                if tan2_theta.is_nan() {
                    0.0
                } else {
                    let cos_phi = spherical_utils::cos_phi(w);
                    let sin_phi = spherical_utils::sin_phi(w);

                    let alpha2 = (cfg.roughness * cos_phi).powi(2)
                        + (cfg.roughness * sin_phi).powi(2);

                    ((1.0 + alpha2 * tan2_theta).max(0.0).sqrt() - 1.0) / 2.0
                }
            }
            Self::Beckmann(cfg) => {
                let tan2_theta = spherical_utils::tan2_theta(w);
                if tan2_theta.is_nan() {
                    0.0
                } else {
                    let a = 1.0 / (cfg.roughness * tan2_theta);

                    if a >= 1.6 {
                        0.0
                    } else {
                        (1.0 - 1.259 * a + 0.396 * a * a)
                            / (3.535 * a + 2.181 * a * a)
                    }
                }
            }
        }
    }

    /// Probability that `wh` got sampled. `wh` and `v` in shading space.
    pub fn sample_normal_pdf(
        &self,
        wh: Normal,
        v: Direction,
    ) -> Float {
        let pdf = match self {
            Self::Beckmann(..) => {
                let cos_theta_wh = spherical_utils::cos_theta(wh);
                self.d(wh) * cos_theta_wh
            }
            Self::Ggx(..) => {
                let wh_dot_v = wh.dot(v);
                let cos_theta_v = spherical_utils::cos_theta(v);

                self.g1(v, wh) * self.d(wh) * wh_dot_v.abs() / cos_theta_v.abs()
            }
        };

        pdf.max(0.0)
    }

    /// Sampling microfacet normals per distribution for importance sampling.
    /// `v` in shading space.
    pub fn sample_normal(&self, v: Direction, rand_sq: Vec2) -> Normal {
        match self {
            Self::Ggx(cfg) => {
                // Heitz 2018 or
                // https://schuttejoe.github.io/post/ggximportancesamplingpart2/

                let roughness = cfg.roughness;
                // Map the GGX ellipsoid to a hemisphere
                let v_stretch = Normal::new(
                    v.x * roughness,
                    v.y * roughness,
                    v.z
                ).normalize();
                if v_stretch.z < 0.0 { -v_stretch } else { v_stretch };

                // ONB basis of the hemisphere configuration
                let u = if v_stretch.z < 0.9999 {
                    v_stretch.cross(Normal::Z).normalize()
                } else {
                    Normal::X
                };
                let vv = v_stretch.cross(u);
                let hemi_basis = Onb::new_from_basis(u, vv, v_stretch);

                // first a point on the unit disk
                let r = rand_sq.x.sqrt();
                let theta = 2.0 * crate::PI * rand_sq.y;
                let x = r * theta.cos();
                // then map it to the projection disk
                let h = (1.0 - x * x).sqrt();
                let lerp = (1.0 + v_stretch.z) / 2.0;
                let y = (1.0 - lerp) * h + lerp * r * theta.sin();

                // compute normal in hemisphere configuration
                let wm = Normal::new(
                    x,
                    y,
                    (1.0 - x*x - y*y).max(0.0).sqrt(),
                );

                // move back to ellipsoid
                let wm = hemi_basis.to_world(wm);
                Normal::new(
                    roughness * wm.x,
                    roughness * wm.y,
                    wm.z.max(1e-6)
                ).normalize()
            }
            Self::Beckmann(cfg) => {
                let roughness2 = cfg.roughness * cfg.roughness;
                let theta = (-roughness2 * (1.0 - rand_sq.y).ln()).sqrt().atan();
                let phi = 2.0 * crate::PI * rand_sq.x;

                // already normalized?
                Normal::new(
                    theta.sin() * phi.cos(),
                    theta.sin() * phi.sin(),
                    theta.cos(),
                ).normalize()
            }
        }
    }
}
