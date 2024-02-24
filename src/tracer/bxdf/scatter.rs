use super::*;

pub fn transmission_f(wo: Direction, eta: Float, mode: Transport) -> Color {
    match mode {
        Transport::Importance => Color::WHITE,
        Transport::Radiance => {
            let inside = wo.z > 0.0;
            if inside {
                Color::splat(1.0 / (eta * eta))
            } else {
                Color::splat(eta * eta)
            }
        }
    }
}

pub fn transmission_sample(wo: Direction, eta: Float) -> Option<Direction> {
    let v = -wo;
    let cos_to = v.z;
    let cos_to_abs = cos_to.abs();
    let inside = cos_to < 0.0;
    let eta_ratio = if inside { eta } else { 1.0 / eta };

    let sin2_to = 1.0 - cos_to * cos_to;
    let sin2_ti = eta_ratio * eta_ratio * sin2_to;

    if sin2_ti >= 1.0 {
        /* total internal reflection */
        // we don't handle it?
        None
    } else {
        let cos_ti = (1.0 - sin2_ti).sqrt();
        let n = if inside { Normal::NEG_Z } else { Normal::Z };
        Some( -v * eta_ratio + (eta_ratio * cos_to_abs - cos_ti) * n )
    }
}

pub fn lambertian_pdf(wi: Direction) -> Float {
    let cos_theta = wi.z;
    if cos_theta > 0.0 {
        cos_theta / crate::PI
    } else {
        0.0
    }
}
