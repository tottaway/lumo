use super::*;

pub fn lambertian_pdf(wi: Direction) -> Float {
    let cos_theta = wi.z;

    if cos_theta > 0.0 {
        cos_theta / crate::PI
    } else {
        0.0
    }
}
