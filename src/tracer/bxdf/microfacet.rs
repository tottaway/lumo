use super::*;

/*
 * MICROFACET TRANSMISSION
 */

pub fn transmission_f(
    wo: Direction,
    wi: Direction,
    mfd: &MfDistribution,
    albedo: Color,
    mode: Transport,
) -> Color {
    let v = -wo;
    let cos_theta_v = v.z;
    let cos_theta_wi = wi.z;

    /* on the same hemisphere */
    if cos_theta_v * cos_theta_wi > 0.0 {
        return Color::BLACK;
    }

    let eta_ratio = if cos_theta_v < 0.0 {
        1.0 / mfd.eta()
    } else {
        mfd.eta()
    };
    let scale = match mode {
        Transport::Radiance => eta_ratio * eta_ratio,
        Transport::Importance => 1.0,
    };
    let wh = (wi * eta_ratio + v).normalize();
    let wh = if wh.z < 0.0 { -wh } else { wh };

    let wh_dot_v = wh.dot(v);
    let wh_dot_wi = wh.dot(wi);


    /* same hemisphere w.r.t. wh */
    if wh_dot_v * wh_dot_wi > 0.0 {
        return Color::BLACK;
    }

    let d = mfd.d(wh, Normal::Z);
    let f = mfd.f(v, wh, albedo);
    let g = mfd.g(v, wi, wh, Normal::Z);

    scale * albedo * d * (Color::WHITE - f) * g
        * (wh_dot_wi * wh_dot_v / (cos_theta_wi * cos_theta_v)).abs()
        / (eta_ratio * wh_dot_wi + wh_dot_v).powi(2)
}

pub fn transmission_sample(
    wo: Direction,
    mfd: &MfDistribution,
    rand_sq: Vec2
) -> Option<Direction> {
    let v = -wo;
    let wh = mfd.sample_normal(v, rand_sq).normalize();
    let inside = v.z < 0.0;
    let eta_ratio = if inside {
        mfd.eta()
    } else {
        1.0 / mfd.eta()
    };

    // here we refract v around wh
    Some( Direction::Z )
}

/*
 * MICROFACET REFLECTION
 */

pub fn reflection_f(
    wo: Direction,
    wi: Direction,
    mfd: &MfDistribution,
    albedo: Color,
) -> Color {
    let v = -wo;
    let cos_theta_v = v.z;
    let cos_theta_wi = wi.z;
    let wh = (wi + v).normalize();

    let d = mfd.d(wh, Normal::Z);
    let f = mfd.f(v, wh, albedo);
    let g = mfd.g(v, wi, wh, Normal::Z);

    d * f * g / (4.0 * cos_theta_v * cos_theta_wi)
}

pub fn reflection_sample(
    wo: Direction,
    mfd: &MfDistribution,
    rand_sq: Vec2
) -> Option<Direction> {
    let v = -wo;
    let wh = self.mfd.sample_normal(v, rand_sq).normalize();
    // reflect v around wh
    let wi = Direction::Z;

    if wi.z <= 0.0 {
        // bad sample, do something else?
        None
    } else {
        Some( wi )
    }

}
