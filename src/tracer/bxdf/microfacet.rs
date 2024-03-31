use super::*;

/* util functions */
mod util {
    use super::*;

    pub fn reflect(v: Direction, no: Normal) -> Option<Direction> {
        let wi = 2.0 * v.project_onto(no) - v;
        if !spherical_utils::same_hemisphere(wi, v) {
            // bad sample, do something else?
            None
        } else {
            Some( wi )
        }
    }

    pub fn refract(eta: Float, v: Direction, no: Normal) -> Option<Direction> {
        /* Snell-Descartes law */
        // possibly flip the orientation if we are "inside"
        let (cos_to, eta_ratio, n) = if no.dot(v) < 0.0 {
            (-no.dot(v), 1.0 / eta, -no)
        } else {
            (no.dot(v), eta, no)
        };
        let sin2_to = 1.0 - cos_to * cos_to;
        let sin2_ti = sin2_to / eta_ratio.powi(2);

        if sin2_ti >= 1.0 {
            /* total internal reflection */
            // we don't do it?
            None
        } else {
            let cos_ti = (1.0 - sin2_ti).max(0.0).sqrt();
            let wi = -v / eta_ratio + (cos_to / eta_ratio - cos_ti) * n;

            if spherical_utils::same_hemisphere(wi, v) {
                None
            } else {
                Some( wi )
            }
        }
    }

    pub fn reflect_coeff(
        wo: Direction,
        wi: Direction,
        mfd: &MfDistribution,
    ) -> Float {
        let v = -wo;
        let cos_theta_v = v.z;
        let cos_theta_wi = wi.z;
        let wh = (wi + v).normalize();

        let d = mfd.d(wh);
        let f = mfd.f(v, wh);
        let g = mfd.g(v, wi, wh);

        // need reflection color, its in the .mtl files somewhere
        d * f * g / (4.0 * cos_theta_v.abs() * cos_theta_wi.abs())
    }
}

/*
 * MICROFACET CONDUCTOR
 */
pub fn conductor_f(
    wo: Direction,
    wi: Direction,
    mfd: &MfDistribution,
    albedo: Color
) -> Color {
    let v = -wo;
    if mfd.is_delta() {
        let f = mfd.f(v, Normal::Z);
        albedo * f / spherical_utils::cos_theta(wi).abs()
    } else {
        albedo * util::reflect_coeff(wo, wi, mfd)
    }
}

pub fn conductor_sample(
    wo: Direction,
    mfd: &MfDistribution,
    rand_sq: Vec2
) -> Option<Direction> {
    let v = -wo;

    if mfd.is_delta() {
        // 2.0 * v.project(Z) - v = 2.0 * (0, 0, v.z) - v = (-v.x, -v.y, v.z)
        Some( Direction::new(-v.x, -v.y, v.z) )
    } else {
        let wh = mfd.sample_normal(v, rand_sq);
        util::reflect(v, wh)
    }
}

pub fn conductor_pdf(
    wo: Direction,
    wi: Direction,
    mfd: &MfDistribution,
) -> Float {
    let v = -wo;

    // check if in same hemisphere or perpendicular to normal
    if !spherical_utils::same_hemisphere(wi, v) {
        return 0.0;
    }

    let wh = (v + wi).normalize();
    let wh = if wh.z < 0.0 { -wh } else { wh };

    if mfd.is_delta() {
        if 1.0 - wh.z < crate::EPSILON { 1.0 } else { 0.0 }
    } else {
        let wh_dot_v = v.dot(wh);
        mfd.sample_normal_pdf(wh, v) / (4.0 * wh_dot_v.abs())
    }
}

/*
 * MICROFACET DIFFUSE
 * Just lambertian
 */
pub fn diffuse_f(
    wo: Direction,
    wi: Direction,
    mfd: &MfDistribution,
    albedo: Color,
) -> Color {
    let v = -wo;
    let wh = (v + wi).normalize();

    /*
    let cos_theta_v = v.z;
    let cos_theta_wi = wi.z;
    let cos_theta_wh = wh.z;
    let disney = mfd.disney_diffuse(cos_theta_v, cos_theta_wi, cos_theta_wh);
    */
    let f = mfd.f(v, wh);

    albedo * ((1.0 - f) / crate::PI + util::reflect_coeff(wo, wi, mfd))
}

/*
 * MICROFACET DIELECTRIC
 */
pub fn dielectric_f(
    wo: Direction,
    wi: Direction,
    mfd: &MfDistribution,
    albedo: Color,
    mode: Transport,
) -> Color {
    let v = -wo;

    let cos_theta_v = spherical_utils::cos_theta(v);
    let cos_theta_wi = spherical_utils::cos_theta(wi);
    let is_reflection = spherical_utils::same_hemisphere(v, wi);
    let v_inside = cos_theta_v < 0.0;

    let eta_ratio = if is_reflection {
        1.0
    } else if v_inside {
        1.0 / mfd.eta()
    } else {
        mfd.eta()
    };

    let wh = if mfd.is_delta() {
        Normal::Z
    } else {
        (wi * eta_ratio + v).normalize()
    };
    let wh = if wh.z < 0.0 { -wh } else { wh };

    let wh_dot_v = wh.dot(v);
    let wh_dot_wi = wh.dot(wi);

    let d = mfd.d(wh);
    let f = mfd.f(v, wh);
    let g = mfd.g(v, wi, wh);

    if is_reflection {
        if mfd.is_delta() {
            albedo * f / cos_theta_wi.abs()
        } else {
            albedo * util::reflect_coeff(wo, wi, mfd)
        }
    } else {
        // scale coefficient if transporting radiance
        let scale = match mode {
            Transport::Radiance => eta_ratio * eta_ratio,
            Transport::Importance => 1.0,
        };

        if mfd.is_delta() {
            albedo * (1.0 - f) / (scale * cos_theta_wi.abs())
        } else {
            albedo * d * (1.0 - f) * g / scale
                * (wh_dot_wi * wh_dot_v / (cos_theta_wi * cos_theta_v)).abs()
                / (eta_ratio * wh_dot_wi + wh_dot_v).powi(2)
        }
    }
}

pub fn dielectric_sample(
    wo: Direction,
    mfd: &MfDistribution,
    rand_sq: Vec2
) -> Option<Direction> {
    let v = -wo;
    let wh = if mfd.is_delta() {
        Normal::Z
    } else {
        mfd.sample_normal(v, rand_sq)
    };

    // importance sample reflection/transmission with fresnel
    let pr = mfd.f(v, wh);
    let pt = 1.0 - pr;

    if rand_utils::rand_float() < pr / (pr + pt) {
        util::reflect(v, wh)
    } else {
        util::refract(mfd.eta(), v, wh)
    }
}

pub fn dielectric_pdf(
    wo: Direction,
    wi: Direction,
    mfd: &MfDistribution,
    swap_dir: bool,
) -> Float {
    let v = -wo;
    let (v, wi) = if swap_dir { (wi, v) } else { (v, wi) };
    let cos_theta_v = spherical_utils::cos_theta(v);
    let cos_theta_wi = spherical_utils::cos_theta(wi);
    let is_reflection = spherical_utils::same_hemisphere(v, wi);
    let v_inside = cos_theta_v < 0.0;

    let eta_ratio = if is_reflection {
        1.0
    } else if v_inside {
        1.0 / mfd.eta()
    } else {
        mfd.eta()
    };

    let wh = (v + wi * eta_ratio).normalize();
    // orient MS normal to same side as geometric normal
    let wh = if wh.z < 0.0 { -wh } else { wh };
    let wh_dot_v = v.dot(wh);
    let wh_dot_wi = wi.dot(wh);

    if wh_dot_v == 0.0 || wh_dot_wi == 0.0 {
        return 0.0;
    }

    // discard backfacing wh
    if wh_dot_v * cos_theta_v < 0.0 || wh_dot_wi * cos_theta_wi < 0.0 {
        return 0.0;
    }

    let pr = mfd.f(v, wh);
    let pt = 1.0 - pr;

    if is_reflection && mfd.is_delta() {
        // reflection with delta
        if 1.0 - wh.z < crate::EPSILON {
            pr / (pr + pt)
        } else {
            0.0
        }
    } else if is_reflection {
        // reflection with rough surface
        mfd.sample_normal_pdf(wh, v) / (4.0 * wh_dot_v.abs())
            * pr / (pr + pt)
    } else if mfd.is_delta() {
        // transmission with delta
        if 1.0 - wh.z < crate::EPSILON {
            pt / (pr + pt)
        } else {
            0.0
        }
    } else {
        // transmission with rough surface
        mfd.sample_normal_pdf(wh, v)
            * wh_dot_wi.abs() / (wh_dot_wi + wh_dot_v / eta_ratio).powi(2)
            * pt / (pr + pt)
    }
}
