use super::*;

/// Mirror with perfect reflection
pub struct Mirror {
    /// Reflected direction
    wi: Direction,
}

impl Mirror {
    pub fn new(ro: &Ray, ho: &Hit) -> Self {
        let v = -ro.dir;
        let no = ho.ns;
        let wi = reflect(v, no);

        Self { xo, wi }
    }
}

impl BSDF for Mirror {
    fn eval(&self, _wo: Direction, _wi: Direction, _mode: Transport) -> Color {
        Color::WHITE
    }

    fn sample_direction(&self, _rand_sq: Vec2) -> Option<Direction> {
        Some( self.wi )
    }

    fn pdf_for(&self, ri: &Ray, swap_dir: bool) -> Float {
        let wi = ri.dir;
        if wi.dot(self.wi) >= 1.0 - crate::EPSILON { 1.0 } else { 0.0 }
    }
}
