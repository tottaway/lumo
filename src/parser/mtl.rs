use super::*;
use crate::tracer::Color;

/// Holds the properties of a microfacet material
#[allow(non_snake_case)]
pub struct MtlConfig {
    /// Diffuse color of the material
    pub Kd: Color,
    /// Specular color of the material
    pub Ks: Color,
    /// Emittance of the material. If not zero vector, then createas a light
    pub Ke: Color,
    /// How much each light channel passes on transmission. Unused.
    pub Tf: Vec3,
    /// Refraction index of the material
    pub Ni: Float,
    /// Roughness of the material
    pub Ns: Float,
    /// Illumination model, see docs.
    /// If 6 or 7 makes transparent, if 5 makes metal, otherwise unused.
    pub illum: usize,
    /// Texture map
    pub map_Kd: Option<Image>,
}

impl Default for MtlConfig {
    fn default() -> Self {
        Self {
            Kd: Color::BLACK,
            Ks: Color::BLACK,
            Ke: Color::BLACK,
            Tf: Vec3::ZERO,
            Ni: 1.5,
            Ns: 0.0,
            illum: 0,
            map_Kd: None,
        }
    }
}

impl MtlConfig {
    pub fn build_material(&self) -> Material {
        if !self.Ke.is_black() {
            Material::Light(Texture::Solid(self.Ke))
        } else {
            let diffuse_color = if let Some(img) = &self.map_Kd {
                Texture::Image(img.clone())
            } else {
                // TODO: separate Ks
                Texture::Solid(self.Kd + self.Ks)
            };
            let _specular_color = Texture::Solid(self.Ks);

            let fresnel_enabled = self.illum == 5 || self.illum == 7;
            let is_transparent = self.illum == 4 || self.illum == 6 || self.illum == 7;
            // blender uses this mapping
            let roughness = 1.0 - self.Ns.min(900.0).sqrt() / 30.0;

            Material::microfacet(
                diffuse_color,
                roughness,
                self.Ni,
                0.0,
                is_transparent,
                fresnel_enabled,
            )
        }
    }
}

pub fn load_file(
    file: File,
    zip_file: Option<Vec<u8>>,
    materials: &mut HashMap<String, MtlConfig>,
) -> Result<()> {
    let reader = BufReader::new(file);

    let mut mtl = MtlConfig::default();
    let mut mtl_name = String::default();

    for line in reader.lines() {
        let line = line?.trim().to_string();
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let tokens: Vec<&str> = line.split_ascii_whitespace().collect();

        match tokens[0] {
            "newmtl" => {
                if !mtl_name.is_empty() {
                    materials.insert(mtl_name, mtl);
                }
                mtl = MtlConfig::default();
                mtl_name = tokens[1].to_string();
            }
            /* diffuse color */
            "Kd" => {
                let kd = parse_vec3(&tokens)?;
                mtl.Kd = Color::from(kd);
            }
            /* texture map */
            "map_Kd" => {
                if let Some(ref zip) = zip_file {
                    let tex_name = tokens[1].replace('\\', "/");
                    let img = super::_img_from_zip(zip.clone(), &tex_name)?;
                    mtl.map_Kd = Some(img);
                }
            }
            /* emission color */
            "Ke" => {
                let ke = parse_vec3(&tokens)?;
                mtl.Ke = Color::from(ke);
            }
            /* specular color */
            "Ks" => {
                let ks = parse_vec3(&tokens)?;
                mtl.Ks = Color::from(ks);
            }
            /* transmission filter */
            "Tf" => {
                let tf = parse_vec3(&tokens)?;
                mtl.Tf = tf;
            }
            /* refraction index */
            "Ni" => {
                let ni = parse_double(tokens[1])?;
                mtl.Ni = ni;
            }
            /* roughness */
            "Ns" => {
                let ns = parse_double(tokens[1])?;
                mtl.Ns = ns;
            }
            /* illumination model */
            "illum" => {
                let illum = parse_double(tokens[1])?;
                mtl.illum = illum as usize;
            }
            _ => (),
        }
    }

    materials.insert(mtl_name, mtl);

    Ok(())
}
