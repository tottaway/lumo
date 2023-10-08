use crate::{ Direction, Normal, Transport, Float, Vec2, rand_utils };
use crate::tracer::{ Color, ray::Ray };

pub enum BxDF {
    Lambertian(Color),
    Transmission(Float),
    Reflection,
}

impl BxDF {
    pub fn eval(&self, wo: Direction, wi: Direction, mode: Transport) -> Color {
        match self {
            Self::Reflection => Color::WHITE,
            Self::Lambertian(albedo) => *albedo / crate::PI,
            Self::Transmission(eta) => {
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
        }
    }

    pub fn sample(&self, wo: Direction, rand_sq: Vec2) -> Option<Direction> {
        match self {
            Self::Reflection => Some( Direction::new(wo.x, wo.y, -wo.z) ),
            Self::Lambertian(_) => {
                Some( rand_utils::square_to_cos_hemisphere(rand_sq) )
            }
            Self::Transmission(eta) => {
                let inside = wo.z > 0.0;
                let eta_ratio = if inside { *eta } else { 1.0 / eta };
                let v = -wo;
                let cos_to = v.z;
                let sin2_to = 1.0 - cos_to * cos_to;
                let sin2_ti = eta_ratio * eta_ratio * sin2_to;

                if sin2_ti > 1.0 {
                    /* total internal reflection */
                    Some( Direction::new(wo.x, wo.y, -wo.z) )
                } else {
                    let cos_ti = (1.0 - sin2_ti).sqrt();
                    let n = if inside { Normal::NEG_Z } else { Normal::Z };
                    Some( -v * eta_ratio + (eta_ratio * cos_to - cos_ti) * n )
                }
            }
        }
    }

    pub fn pdf(&self, wo: Direction, wi: Direction, swap_dir: bool) -> Float {
        match self {
            Self::Reflection | Self::Transmission(_) => 1.0,
            Self::Lambertian(_) => {
                let cos_theta = wi.z;
                if cos_theta > 0.0 {
                    cos_theta / crate::PI
                } else {
                    0.0
                }
            }
        }
    }
}
