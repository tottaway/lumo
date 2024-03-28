use crate::{ Direction, Float };

#[cfg(test)]
mod spherical_utils_tests;

pub fn phi(w: Direction) -> Float {
    let phi = w.y.atan2(w.x);
    if phi < 0.0 { phi + 2.0 * crate::PI } else { phi }
}

pub fn cos_phi(w: Direction) -> Float {
    let sin_theta = sin_theta(w);
    if sin_theta == 0.0 { 1.0 } else { (w.x / sin_theta).clamp(-1.0, 1.0) }
}

pub fn sin_phi(w: Direction) -> Float {
    let sin_theta = sin_theta(w);
    if sin_theta == 0.0 { 0.0 } else { (w.y / sin_theta).clamp(-1.0, 1.0) }
}

pub fn theta(w: Direction) -> Float {
    w.z.clamp(-1.0, 1.0).acos()
}

pub fn cos_theta(w: Direction) -> Float {
    w.z
}

pub fn sin_theta(w: Direction) -> Float {
    sin2_theta(w).sqrt()
}

pub fn tan_theta(w: Direction) -> Float {
    let cos = cos_theta(w);

    if cos < crate::EPSILON {
        Float::NAN
    } else {
        sin_theta(w) / cos_theta(w)
    }
}

pub fn cos2_theta(w: Direction) -> Float {
    w.z * w.z
}

pub fn sin2_theta(w: Direction) -> Float {
    (1.0 - cos2_theta(w)).max(0.0)
}

pub fn tan2_theta(w: Direction) -> Float {
    let cos2 = cos2_theta(w);

    if cos2 < crate::EPSILON {
        Float::NAN
    } else {
        sin2_theta(w) / cos2_theta(w)
    }
}
