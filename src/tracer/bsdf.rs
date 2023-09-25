use crate::{ Direction, Point, Transport, Float, Vec2 };
use crate::tracer::{ Color, ray::Ray };

pub trait BSDF {
    /// Evaluate BSDF for incoming `wo` and outgoing `wi` with `mode`
    fn eval(&self, wo: Direction, wi: Direction, mode: Transport) -> Color;
    /// Sample direction given a random value in unit square
    fn sample_direction(&self, rand_sq: Vec2) -> Option<Direction>;
    /// Probability that `ri` got sampled. `swap_dir` if inverse transport
    fn pdf_for(&self, ri: &Ray, swap_dir: bool) -> Float;
}


/// Reflect around normal
///
/// # Arguments
/// * `v` - Normalized? direction from reflection point to viewer
/// * `no` - Surface normal
fn reflect(v: Direction, no: Normal) -> Direction {
    2.0 * v.project_onto(no) - v
}

/// Refract direction with Snell-Descartes law.
///
/// # Arguments
/// * `eta_ratio` - Ratio of refraction indices. `from / to`
/// * `v` - Normalized direction from refraction point to viewer
/// * `no` - Surface normal, pointing to same hemisphere as `v`
fn refract(eta_ratio: Float, v: Direction, no: Normal) -> Direction {
    /* Snell-Descartes law */
    let cos_to = no.dot(v);
    let sin2_to = 1.0 - cos_to * cos_to;
    let sin2_ti = eta_ratio * eta_ratio * sin2_to;

    /* total internal reflection */
    if sin2_ti > 1.0 {
        reflect(v, no)
    } else {
        let cos_ti = (1.0 - sin2_ti).sqrt();

        -v * eta_ratio + (eta_ratio * cos_to - cos_ti) * no
    }
}
