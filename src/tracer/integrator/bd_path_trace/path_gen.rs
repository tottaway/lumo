use super::*;

/// Generates a ray path starting from the camera
pub fn camera_path<'a>(scene: &'a Scene, camera: &'a Camera, r: Ray) -> Vec<Vertex<'a>> {
    let gathered = Color::WHITE;
    let root = Vertex::camera(r.origin, gathered);
    let wi = r.dir;
    let pdf_fwd = camera.pdf_wi(wi);

    walk(scene, r, root, gathered, pdf_fwd, Transport::Radiance)
}

/// Generates a ray path strating from a light
pub fn light_path(scene: &Scene) -> Vec<Vertex> {
    let light = scene.uniform_random_light();
    let pdf_light = 1.0 / scene.num_lights() as Float;
    let (ro, ho) = light.sample_leaving(
        rand_utils::unit_square(),
        rand_utils::unit_square()
    );
    let ng = ho.ng;
    let ns = ho.ns;
    let (pdf_origin, pdf_dir) = light.sample_leaving_pdf(&ro, ng);
    let emit = ho.material.emit(&ho);
    let root = Vertex::light(ho, light, emit, pdf_origin * pdf_light);

    let gathered = emit * ns.dot(ro.dir).abs()
        / (pdf_light * pdf_origin * pdf_dir);

    walk(scene, ro, root, gathered, pdf_dir, Transport::Importance)
}

/// Ray that randomly scatters around from the given root vertex
fn walk<'a>(
    scene: &'a Scene,
    mut ro: Ray,
    root: Vertex<'a>,
    mut gathered: Color,
    pdf_dir: Float,
    mode: Transport,
) -> Vec<Vertex<'a>> {
    let mut depth = 0;
    let mut vertices = vec![root];
    let mut pdf_fwd = pdf_dir;

    while let Some(ho) = scene.hit(&ro) {
        let material = ho.material;
        gathered *= scene.transmittance(ho.t);

        let prev = depth;
        let curr = depth + 1;
        vertices.push(Vertex::surface(
            ho,
            gathered,
            pdf_fwd,
            &vertices[prev],
        ));
        let ho = &vertices[curr].h;
        let wo = ro.dir;
        match material.bsdf_sample(wo, ho, rand_utils::unit_square()) {
            None => {
                // we hit a light. if tracing from a light, discard latest vertex
                if matches!(mode, Transport::Importance) {
                    vertices.pop();
                }
                break;
            }
            Some(wi) => {
                let xo = ho.p;
                let wo = ro.dir;
                let ng = ho.ng;

                let ns = ho.ns;
                let ri = ho.generate_ray(wi);
                // normalized
                let wi = ri.dir;

                pdf_fwd = material.bsdf_pdf(wo, wi, ho, false);

                if pdf_fwd <= 0.0 {
                    break;
                }

                let shading_cosine = match mode {
                    Transport::Radiance => material.shading_cosine(wi, ns),
                    Transport::Importance => {
                        if ho.is_medium() {
                            1.0
                        } else {
                            let xp = vertices[prev].h.p;
                            let v = (xp - xo).normalize();
                            wi.dot(ng).abs()
                                * material.shading_cosine(v, ns)
                                / v.dot(ng).abs()
                        }
                    }
                };

                let bsdf = material.bsdf_f(wo, wi, mode, ho);
                let bsdf = if ho.is_medium() {
                    bsdf * pdf_fwd
                } else {
                    bsdf
                };

                gathered *= bsdf * shading_cosine / pdf_fwd;

                vertices[prev].pdf_bck = if material.is_delta() || !vertices[prev].is_surface() {
                    0.0
                } else {
                    let pdf_bck = material.bsdf_pdf(wo, wi, ho, true);
                    vertices[curr].solid_angle_to_area(pdf_bck, &vertices[prev])
                };

                if material.is_delta() {
                    pdf_fwd = 0.0;
                }

                // russian roulette
                if depth > 3 {
                    let luminance = gathered.luminance();
                    let rr_prob = (1.0 - luminance).max(0.05);
                    if rand_utils::rand_float() < rr_prob {
                        break;
                    }

                    // TODO (9)
                    //gathered /= 1.0 - rr_prob;
                }

                depth += 1;
                ro = ri;
            }
        }
    }

    vertices
}
