use super::*;

const NUM_BINS: usize = 5;
const NUM_SAMPLES: usize = 10_000_000;

#[test]
fn lambertian_sampling() {
    let bxdf = BxDF::Lambertian;

    assert_bins(do_sampling(bxdf));
}

#[test]
fn conductor_sampling_75() {
    let mfd = MfDistribution::new(0.75, 1.0, 0.0, true);
    let bxdf = BxDF::MfConductor(mfd);

    assert_bins(do_sampling(bxdf));
}

#[test]
fn conductor_sampling_50() {
    let mfd = MfDistribution::new(0.5, 1.0, 0.0, true);
    let bxdf = BxDF::MfConductor(mfd);

    assert_bins(do_sampling(bxdf));
}

#[test]
fn conductor_sampling_25() {
    let mfd = MfDistribution::new(0.25, 1.0, 0.0, true);
    let bxdf = BxDF::MfConductor(mfd);

    assert_bins(do_sampling(bxdf));
}

#[test]
fn conductor_sampling_10() {
    let mfd = MfDistribution::new(0.10, 1.0, 0.0, true);
    let bxdf = BxDF::MfConductor(mfd);

    assert_bins(do_sampling(bxdf));
}

fn print_bins(bins: &Vec<Vec<Float>>) {
    println!("bin values:");
    for i in 0..NUM_BINS {
        for j in 0..NUM_BINS {
            print!("{:.3} ", bins[i][j]);
        }
        println!("");
    }
}

fn assert_bins(bins: Vec<Vec<Float>>) {
    // print bins first for debugging
    print_bins(&bins);
    let mut sum = 0.0;
    println!("bin deltas:");
    for i in 0..NUM_BINS {
        for j in 0..NUM_BINS {
            sum += bins[i][j];
            let delta = (2.0 * crate::PI - bins[i][j]).abs();
            print!("{:.3} ", delta);
            // rather lax here...
            assert!(delta < NUM_BINS as Float * 1e-1);
        }
        println!("");
    }

    sum /= (NUM_BINS * NUM_BINS) as Float;
    let delta = (2.0 * crate::PI - sum).abs();
    println!("delta: {}", delta);
    assert!(delta < 1e-2);
}

/* ripped from pbrt. take v = normal and sample NUM_SAMPLES times.
 * if sample is ok add 1 / pdf to a *bin*. each bin should converge to 2*PI.
 * split bins using spherical coordinates of the sampled vector.
 */
fn do_sampling(bxdf: BxDF) -> Vec<Vec<Float>> {
    let mut bins: Vec<Vec<Float>> = vec!();
    let mut num_failed: usize = 0;

    for i in 0..NUM_BINS {
        bins.push(vec!());
        for _ in 0..NUM_BINS {
            bins[i].push(0.0);
        }
    }
    // let v face directly the normal
    let wo = Direction::NEG_Z;
    for _ in 0..NUM_SAMPLES {
        let wi = bxdf.sample(wo, rand_utils::unit_square());
        match wi {
            None => num_failed += 1,
            Some(wi) => {
                let pdf = bxdf.pdf(wo, wi, false);
                if pdf == 0.0 {
                    num_failed += 1;
                    continue;
                }
                let wi_phi = (wi.x.atan2(wi.y) + crate::PI) / (2.0 * crate::PI);
                let phi_idx = (wi_phi * NUM_BINS as Float) as usize;
                let wi_cos_theta = wi.z;
                let theta_idx = (wi_cos_theta * NUM_BINS as Float) as usize;
                // can be out of bounds in rare cases. fix if it happens.
                bins[theta_idx][phi_idx] += 1.0 / pdf;
            }
        }
    }
    let good_samples = NUM_SAMPLES - num_failed;

    // scale bins properly
    for i in 0..NUM_BINS {
        for j in 0..NUM_BINS {
            bins[i][j] *= (NUM_BINS * NUM_BINS) as Float / good_samples as Float;
        }
    }

    bins
}
