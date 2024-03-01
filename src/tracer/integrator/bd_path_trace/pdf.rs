use super::*;
/// Randomly samples a direction towards a point on the object that is visible
pub struct ObjectPdf<'a> {
    /// Object to do sampling from
    object: &'a dyn Sampleable,
    /// Point from where the object should be visible
    xo: Point,
}

impl<'a> ObjectPdf<'a> {
    pub fn new(object: &'a dyn Sampleable, xo: Point) -> Self {
        Self { object, xo }
    }

    pub fn sample_direction(&self, rand_sq: Vec2) -> Option<Direction> {
        Some( self.object.sample_towards(self.xo, rand_sq) )
    }

    pub fn value_for(&self, ri: &Ray, _swap_dir: bool) -> Float {
        let (p, hi) = self.object.sample_towards_pdf(ri);
        if let Some(hi) = hi {
            // convert area measure to solid angle measure
            // other fields of hit might be in local instance coordinates
            let xi = hi.p;
            let ni = hi.ng;
            let wi = ri.dir;
            p * self.xo.distance_squared(xi) / ni.dot(wi).abs()
        } else {
            0.0
        }
    }
}
