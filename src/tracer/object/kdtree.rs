use super::*;
use crate::cli::TracerCli;
use clap::Parser;
use std::time::Instant;

/// Triangle mesh constructed as a kD-tree
pub type Mesh = KdTree<Triangle>;

#[cfg(test)]
mod kdtree_tests;

/// A k dimensional tree used to accelerate ray intersection calculations.
/// Implements a binary tree that splits a large mesh of objects to smaller
/// subobjects.
pub struct KdTree<T> {
    objects: Vec<T>,
    boundary: AaBoundingBox,
    root: Box<KdNode>,
    material: Material,
}

/// Implementation of a SAH based kd-tree
/// References:
/// [ekzhang/rpt](https://github.com/ekzhang/rpt/blob/master/src/kdtree.rs),
/// [fogleman/pt](https://github.com/fogleman/pt/blob/master/pt/tree.go),
/// [Article by Amsallem](https://www.flomonster.fr/articles/kdtree.html)
impl<T: Bounded> KdTree<T> {
    /// Constructs a kD-tree of the given objects with the given material.
    /// Should each object have their own material instead?
    pub fn new(objects: Vec<T>, material: Material) -> Self {
        let start = Instant::now();
        if objects.len() > 10_000 {
            println!("Creating kd-tree of {} triangles", objects.len());
        }

        let indices = (0..objects.len()).collect();
        let bounds: Vec<AaBoundingBox> = objects.iter().map(|obj| obj.bounding_box()).collect();
        let boundary = bounds
            .iter()
            .fold(AaBoundingBox::default(), |b1, b2| b1.merge(b2));

        let threads = match TracerCli::try_parse() {
            Ok(cli_args) => cli_args.threads.unwrap_or(0),
            Err(_) => 0,
        };
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        let root = pool.install(|| KdNode::construct(&bounds, &boundary, indices));

        if objects.len() > 10_000 {
            println!("Created kd-tree in {:#?}", start.elapsed());
        }

        Self {
            root,
            objects,
            boundary,
            material,
        }
    }

    /// Returns self uniformly scaled as an instance with largest dimension
    /// of bounding box scaled to 1.0
    pub fn to_unit_size(self) -> Box<Instance<Self>> {
        let AaBoundingBox { ax_min, ax_max } = self.bounding_box();

        let bb_dim = ax_max - ax_min;
        let s = 1.0 / bb_dim.max_element();
        self.scale(s, s, s)
    }

    fn hit_subtree(
        &self,
        node: &KdNode,
        r: &Ray,
        t_min: Float,
        t_max: Float,
        aabb: &AaBoundingBox,
    ) -> Option<Hit> {
        // extract split info or check for hit at leaf node
        let (axis, point, mut node_first, mut node_second) = match node {
            KdNode::Split(axis, point, left, right) => (*axis, *point, left, right),
            KdNode::Leaf(indices) => {
                let mut tt = t_max;
                let mut h = None;
                for idx in indices {
                    h = self.objects[*idx].hit(r, t_min, tt).or(h);
                    tt = h.as_ref().map_or(tt, |hit| hit.t);
                }
                return h.map(|mut h| {
                    h.material = &self.material;
                    h
                });
            }
        };

        let t_split = (point - r.o(axis)) / r.d(axis);

        let (mut aabb_first, mut aabb_second) = aabb.split(axis, point);

        let left_first = r.o(axis) < point || (r.o(axis) == point && r.d(axis) <= 0.0);
        // intersect first the AABB that we reach first
        if !left_first {
            std::mem::swap(&mut aabb_first, &mut aabb_second);
            std::mem::swap(&mut node_first, &mut node_second);
        }

        let (t_start, t_end) = aabb.intersect(r);
        let (t_start, t_end) = (t_start.max(t_min), t_end.min(t_max));

        // PBR Figure 4.19 (a). we hit only the first aabb.
        if t_split > t_end || t_split <= 0.0 {
            self.hit_subtree(node_first, r, t_min, t_end, &aabb_first)
        // PBR Figure 4.19 (b). we hit only the second aabb.
        } else if t_split < t_start {
            self.hit_subtree(node_second, r, t_min, t_end, &aabb_second)
        } else {
            match self.hit_subtree(node_first, r, t_start, t_end, &aabb_first) {
                None => self.hit_subtree(node_second, r, t_min, t_end, &aabb_second),
                Some(h1) => {
                    /* if we hit something in the first AABB before the split,
                     * there is no need to process the other subtree. */
                    if h1.t < t_split {
                        Some(h1)
                    } else {
                        self.hit_subtree(node_second, r, t_min, h1.t, &aabb_second)
                            .or(Some(h1))
                    }
                }
            }
        }
    }
}

impl<T: Bounded> Bounded for KdTree<T> {
    fn bounding_box(&self) -> AaBoundingBox {
        self.boundary
    }
}

impl<T: Bounded> Object for KdTree<T> {
    fn hit(&self, r: &Ray, t_min: Float, t_max: Float) -> Option<Hit> {
        let (t_start, t_end) = self.boundary.intersect(r);
        let (t_start, t_end) = (t_start.max(t_min), t_end.min(t_max));
        // box missed / is behind
        if t_start > t_end {
            None
        } else {
            self.hit_subtree(&self.root, r, t_min, t_end, &self.boundary)
        }
    }
}

impl<T: Sampleable + Bounded> Sampleable for KdTree<T> {
    fn area(&self) -> Float {
        // maybe sloooow for big ones
        self.objects.iter().fold(0.0, |sum, obj| sum + obj.area())
    }

    fn sample_on(&self, rand_sq: Vec2) -> Hit {
        let n = rand_utils::rand_float() * self.objects.len() as Float;
        let mut ho = self.objects[n.floor() as usize].sample_on(rand_sq);
        ho.material = &self.material;
        ho
    }
}

const COST_TRAVERSE: Float = 15.0;
const COST_INTERSECT: Float = 20.0;
const EMPTY_BONUS: Float = 0.2;

/// A node in the kD-tree. Can be either a plane split or a leaf node.
pub enum KdNode {
    /// X-split, axis (x = 0, y = 1, z = 2), split point and child nodes
    Split(Axis, Float, Box<KdNode>, Box<KdNode>),
    /// Stores indices to the object vector in the kD-tree
    Leaf(Vec<usize>),
}

impl KdNode {
    /// Computes the cost for split along `axis` at `point`.
    fn cost(
        boundary: &AaBoundingBox,
        axis: Axis,
        point: Float,
        num_left: usize,
        num_right: usize,
    ) -> Float {
        if !boundary.cuts(axis, point) {
            crate::INF
        } else {
            let (left, right) = boundary.split(axis, point);

            let cost = COST_TRAVERSE
                + COST_INTERSECT
                    * (num_left as Float * left.area() / boundary.area()
                        + num_right as Float * right.area() / boundary.area());

            if num_left == 0 || num_right == 0 {
                (1.0 - EMPTY_BONUS) * cost
            } else {
                cost
            }
        }
    }

    /// Finds the best split according to SAH.
    fn find_best_split(
        aabbs: &Vec<&AaBoundingBox>,
        boundary: &AaBoundingBox,
    ) -> (Axis, Float, Float) {
        let mut best_cost = crate::INF;
        let mut best_point = crate::INF;
        let mut best_axis = Axis::X;

        for axis in [Axis::X, Axis::Y, Axis::Z] {
            let mut mins: Vec<Float> = Vec::with_capacity(aabbs.len());
            let mut maxs: Vec<Float> = Vec::with_capacity(aabbs.len());

            aabbs.iter().for_each(|aabb| {
                mins.push(aabb.min(axis));
                maxs.push(aabb.max(axis));
            });

            mins.sort_by(|a: &Float, b: &Float| a.partial_cmp(b).unwrap());
            maxs.sort_by(|a: &Float, b: &Float| a.partial_cmp(b).unwrap());

            let mut num_left = 0;
            let mut num_right = aabbs.len();

            // iterator for mins
            let mut min_idx = 0;
            // iterator for maxs
            let mut max_idx = 0;

            // add infinity to end as "null"
            mins.push(crate::INF);
            maxs.push(crate::INF);

            // do quasi merge
            while mins[min_idx] < crate::INF || maxs[max_idx] < crate::INF {
                let is_min = mins[min_idx] <= maxs[max_idx];
                let point = mins[min_idx].min(maxs[max_idx]);

                // update objects on right before cost..
                if !is_min {
                    max_idx += 1;
                    num_right -= 1;
                }

                let cost = Self::cost(boundary, axis, point, num_left, num_right);

                if cost < best_cost {
                    best_cost = cost;
                    best_axis = axis;
                    best_point = point;
                }

                // ..and objects on left after cost
                if is_min {
                    min_idx += 1;
                    num_left += 1;
                }
            }
        }

        (best_axis, best_point, best_cost)
    }

    /// Partitions indices to left and right parts along `axis` at `point`
    fn partition(
        aabbs: &[&AaBoundingBox],
        indices: Vec<usize>,
        axis: Axis,
        point: Float,
    ) -> (Vec<usize>, Vec<usize>) {
        let mut left: Vec<usize> = Vec::with_capacity(aabbs.len());
        let mut right: Vec<usize> = Vec::with_capacity(aabbs.len());
        aabbs.iter().zip(indices).for_each(|(aabb, idx)| {
            if aabb.min(axis) < point {
                left.push(idx);
            }
            // are we missing one, since both check strict?
            if aabb.max(axis) > point {
                right.push(idx);
            }
        });
        (left, right)
    }

    /// Constructs nodes of the kD-tree recursively with SAH.
    pub fn construct(
        bounds: &[AaBoundingBox],
        boundary: &AaBoundingBox,
        indices: Vec<usize>,
    ) -> Box<Self> {
        // filter relevant AABBs
        let aabbs: Vec<&AaBoundingBox> = indices.iter().map(|idx| &bounds[*idx]).collect();
        let (axis, point, cost) = Self::find_best_split(&aabbs, boundary);
        // cut not worth it, make a leaf
        if cost > COST_INTERSECT * indices.len() as Float {
            Box::new(Self::Leaf(indices))
        } else {
            let (left_idx, right_idx) = Self::partition(&aabbs, indices, axis, point);
            let (left_bound, right_bound) = boundary.split(axis, point);
            let (left, right) = rayon::join(
                || Self::construct(bounds, &left_bound, left_idx),
                || Self::construct(bounds, &right_bound, right_idx),
            );
            Box::new(Self::Split(axis, point, left, right))
        }
    }
}
