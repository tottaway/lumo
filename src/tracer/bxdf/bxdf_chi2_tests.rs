use super::*;
use statrs::distribution::{ ChiSquared, ContinuousCDF };

const NUM_SAMPLES: usize = 1_000_000;
const THETA_BINS: usize = 64;
const PHI_BINS: usize = 2 * THETA_BINS;
const CHI2_RUNS: usize = 5;
const CHI2_SLEVEL: Float = 0.0001;

fn mfd(roughness: Float) -> MfDistribution {
    MfDistribution::new(roughness, 1.0, 0.0, true)
}

#[test]
fn lambertian_chi2() {
    let bxdf = BxDF::Lambertian;
    for _ in 0..CHI2_RUNS {
        let wo = Direction::NEG_Z;
        assert!(chi2_pass(wo, &bxdf));
    }
    assert!(false);
}

#[test]
fn conductor50_chi2() {
    let mfd = mfd(0.5);
    let bxdf = BxDF::MfConductor(mfd);
    for _ in 0..CHI2_RUNS {
        let wo = Direction::NEG_Z;
        assert!(chi2_pass(wo, &bxdf));
    }
    assert!(false);
}

fn chi2_pass(wo: Direction, bxdf: &BxDF) -> bool {
    let actual_freq = sample_frequencies(wo, &bxdf);
    let expected_freq = compute_frequencies(wo, &bxdf);

    // degrees of freedom
    let mut dof = 0;
    // test statistic = sum_i (O_i - D_i)^2 / D_i
    let mut stat = 0.0;

    for phi_bin in 0..PHI_BINS {
        for theta_bin in 0..THETA_BINS {
            let idx = phi_bin + theta_bin * PHI_BINS;
            let delta = actual_freq[idx] as Float - expected_freq[idx] as Float;
            stat += (delta * delta) / expected_freq[idx] as Float;
            dof += 1;
        }
    }

    // blaa blaa i dont know
    dof -= 1;

    // todo: code the CDF, should be simple but boring
    let chi2 = ChiSquared::new(dof as Float).unwrap();
    let pval = 1.0 - chi2.cdf(stat);

    let alpha = 1.0 - CHI2_SLEVEL;
    println!("{} {}", pval, stat);
    pval < alpha
}

fn sample_frequencies(wo: Direction, bxdf: &BxDF) -> [usize; THETA_BINS*PHI_BINS] {
    let mut samples = [0; THETA_BINS*PHI_BINS];

    for _ in 0..NUM_SAMPLES {
        match bxdf.sample(wo, rand_utils::unit_square()) {
            None => (),
            Some(wi) => {
                let theta = wi.z.acos();
                let phi = wi.y.atan2(wi.x);
                let phi = if phi < 0.0 { phi + 2.0 * crate::PI } else { phi };

                let theta_bin = theta * THETA_BINS as Float / crate::PI;
                let theta_bin = (theta_bin as usize).max(0).min(THETA_BINS - 1);
                let phi_bin = phi * PHI_BINS as Float / (2.0 * crate::PI);
                let phi_bin = (phi_bin as usize).max(0).min(PHI_BINS - 1);

                samples[phi_bin + theta_bin * PHI_BINS] += 1;
            }
        }
    }

    return samples;
}

fn compute_frequencies(wo: Direction, bxdf: &BxDF) -> [usize; THETA_BINS*PHI_BINS] {
    let mut samples = [0; THETA_BINS*PHI_BINS];

    let theta_factor = crate::PI / THETA_BINS as Float;
    let phi_factor = (2.0 * crate::PI) / PHI_BINS as Float;

    for theta_bin in 0..THETA_BINS {
        let theta0 = theta_bin as Float * theta_factor;
        let theta1 = theta0 + theta_factor;
        for phi_bin in 0..PHI_BINS {
            let phi0 = phi_bin as Float * phi_factor;
            let phi1 = phi0 + phi_factor;
            let f = |theta: Float, phi: Float| {
                let wi = Direction::new(
                    theta.sin() * phi.cos(),
                    theta.sin() * phi.sin(),
                    theta.sin(),
                );
                bxdf.pdf(wo, wi, false)
            };
            let integral = simpson2d(f, theta0, phi0, theta1, phi1);
            samples[phi_bin + theta_bin * PHI_BINS] = (integral * NUM_SAMPLES as Float) as usize;
        }
    }

    return samples;
}

fn simpson2d<F>(f: F, x0: Float, y0: Float, x1: Float, y1: Float) -> Float where
    F: Fn(Float, Float) -> Float {
    let g = |y: Float| {
        let h = |x: Float| f(x, y);
        simpson(h, x0, x1)
    };
    simpson(g, y0, y1)
}

fn simpson<F>(f: F, x0: Float, x1: Float) -> Float where
    F: Fn(Float) -> Float {
    (f(x0) + 4.0 * f((x0 + x1) / 2.0) + f(x1)) * (x1 - x0) / 6.0
}
