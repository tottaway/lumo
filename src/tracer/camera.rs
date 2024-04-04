use crate::{
    Point, Direction, Float, Vec2, Transform, Normal,
    Mat3, Mat4, Vec4, Vec3, rand_utils, spherical_utils
};
use glam::IVec2;
use crate::tracer::{
    film::FilmSample, ray::Ray, Color
};

/// Common configuration for cameras
pub struct CameraConfig {
    /// Camera position in world space
    pub origin: Point,
    /// Image resolution
    pub resolution: IVec2,
    /// Focal length i.e. distance to focal point behind camera
    pub focal_length: Float,
    /// Radius of the camera lens
    pub lens_radius: Float,
    /// Screen space to camera space transformation
    pub screen_to_camera: Mat4,
    /// Raster space to screen space transformation
    pub raster_to_screen: Transform,
    /// Camera space to world space transformation
    camera_to_world: Transform,
}

impl CameraConfig {
    /// Creates a new config with the given arguments
    pub fn new(
        origin: Point,
        towards: Point,
        up: Direction,
        lens_radius: Float,
        focal_length: Float,
        resolution: (i32, i32),
        screen_to_camera: Mat4,
    ) -> Self {
        assert!(resolution.0 > 0 && resolution.1 > 0);
        assert!(lens_radius >= 0.0);
        assert!(origin != towards);
        assert!(up.length() != 0.0);

        // x = right, y = down, z = towards
        let forward = (towards - origin).normalize();
        let right = forward.cross(up).normalize();
        let down = forward.cross(right);
        let camera_to_world = Transform::from_mat3_translation(
            Mat3::from_cols(right, down, forward).transpose(),
            origin,
        ).inverse();

        let (width, height) = resolution;

        let screen_min = Vec2::new(-1.0, -1.0);
        let screen_max = Vec2::new(1.0, 1.0);
        let screen_delta = screen_max - screen_min;
        let screen_to_raster =
            // ndc_to_raster
            Transform::from_scale(Vec3::new(width as Float, -height as Float, 1.0))
            // screen_to_ndc
            * Transform::from_scale(Vec3::new(1.0 / screen_delta.x, 1.0 / screen_delta.y, 1.0))
            * Transform::from_translation(Vec3::new(-screen_min.x, -screen_max.y, 0.0));
        let raster_to_screen = screen_to_raster.inverse();

        Self {
            lens_radius,
            focal_length,
            origin,
            screen_to_camera,
            camera_to_world,
            raster_to_screen,
            resolution: IVec2::new(width, height),
        }
    }

    pub fn point_to_local(&self, xo: Point) -> Point {
        self.camera_to_world.inverse().transform_point3(xo)
    }

    pub fn point_to_world(&self, xo_local: Point) -> Point {
        self.camera_to_world.transform_point3(xo_local)
    }

    pub fn direction_to_local(&self, wi: Direction) -> Direction {
        self.camera_to_world.inverse().transform_vector3(wi)
    }

    pub fn direction_to_world(&self, wi_local: Direction) -> Direction {
        self.camera_to_world.transform_vector3(wi_local)
    }

    pub fn normal_to_local(&self, no: Normal) -> Normal {
        let m = self.camera_to_world.matrix3.transpose();
        no.x * m.x_axis + no.y * m.y_axis + no.z * m.z_axis
    }

    pub fn normal_to_world(&self, no: Normal) -> Normal {
        let m = self.camera_to_world.matrix3.inverse().transpose();
        no.x * m.x_axis + no.y * m.y_axis + no.z * m.z_axis
    }

    pub fn raster_to_camera(&self, raster_xy: Vec2) -> Point {
        let raster_xyz = raster_xy.extend(0.0);
        let screen_xyz = self.raster_to_screen.transform_point3(raster_xyz);
        (self.screen_to_camera * screen_xyz.extend(1.0)).truncate()
    }

    pub fn camera_to_raster(&self, xo_local: Point) -> Vec2 {
        let screen_xyz = self.screen_to_camera.inverse().project_point3(xo_local);
        let raster_xyz = self.raster_to_screen.inverse().transform_point3(screen_xyz);
        raster_xyz.truncate()
    }
}

/// Camera abstraction
pub enum Camera {
    /// Perspective camera with configurable vertical field-of-view
    Perspective(CameraConfig),
    /// Orthographic camera that preserves angles with configurable image plane scale
    Orthographic(CameraConfig, Float),
}

impl Camera {
    /// Orthographic camera that preserves angles. All rays are cast in the same
    /// direction but from a plane instead of a single point
    ///
    /// # Arguments
    /// * `origin` - Camera position in world space
    /// * `towards` - Point in world space the camera is looking at
    /// * `up` - Up direction of the camera
    /// * `image_plane_scale` - Scale of the plane rays are cast from
    /// * `lens_radius` - Radius of the lens for depth of field. Bigger means more profound effect
    /// * `focal_length` - Distance to the plane of focus for depth of field
    /// * `width` - Width of the rendered image
    /// * `height` - Height of the rendered image
    #[allow(clippy::too_many_arguments)]
    pub fn orthographic(
        origin: Point,
        towards: Point,
        up: Direction,
        image_plane_scale: Float,
        lens_radius: Float,
        focal_length: Float,
        width: i32,
        height: i32,
    ) -> Self {
        assert!(image_plane_scale > 0.0);

        let near = 0.0;
        let far = 1.0;
        let camera_to_screen = Mat4::from_scale(Vec3::new(1.0, 1.0, 1.0 / (far - near)))
            * Mat4::from_translation(Vec3::new(0.0, 0.0, -near));

        Self::Orthographic(
            CameraConfig::new(
                origin,
                towards,
                up,
                lens_radius,
                focal_length,
                (width, height),
                camera_to_screen.inverse(),
            ),
            image_plane_scale,
        )
    }

    /// Perspective camera where sense of depth is more profound. Rays are cast
    /// from a single point towards points on the image plane.
    ///
    /// # Arguments
    /// * `origin` - Camera position in world space
    /// * `towards` - Point in world space the camera is looking at
    /// * `up` - Up direction of the camera
    /// * `vfov` - Vertical field of view of the camera
    /// * `lens_radius` - Radius of the lens for depth of field. Bigger means more profound effect
    /// * `focal_length` - Distance to the plane of focus for depth of field
    /// * `width` - Width of the rendered image
    /// * `height` - Height of the rendered image
    #[allow(clippy::too_many_arguments)]
    pub fn perspective(
        origin: Point,
        towards: Point,
        up: Direction,
        vfov: Float,
        lens_radius: Float,
        focal_length: Float,
        width: i32,
        height: i32,
    ) -> Self {
        assert!(vfov > 0.0 && vfov < 180.0);

        let near = 1e-2;
        let far = 1e3;
        let a = far / (far - near);
        let b = -far * near / (far - near);
        let projection = Mat4::from_cols(
            Vec4::new(1.0, 0.0, 0.0, 0.0),
            Vec4::new(0.0, 1.0, 0.0, 0.0),
            Vec4::new(0.0, 0.0, a,   b),
            Vec4::new(0.0, 0.0, 1.0, 0.0),
        ).transpose();
        let tan_vfov_inv = 1.0 / (vfov.to_radians() / 2.0);
        let scale = Mat4::from_scale(Vec3::new(tan_vfov_inv, tan_vfov_inv, 1.0));
        let camera_to_screen = scale * projection;

        Self::Perspective(
            CameraConfig::new(
                origin,
                towards,
                up,
                lens_radius,
                focal_length,
                (width, height),
                camera_to_screen.inverse(),
            )
        )
    }

    /// The "default" camera. Perspective camera at world space origin
    /// pointing towards `-z` with `y` as up and vfov at 90° with no DOF
    pub fn default(width: i32, height: i32) -> Self {
        Self::perspective(
            Point::ZERO,
            Point::NEG_Z,
            Direction::Y,
            90.0,
            0.0,
            0.0,
            width,
            height,
        )
    }

    fn get_cfg(&self) -> &CameraConfig {
        match self {
            Self::Orthographic(cfg, _) | Self::Perspective(cfg) => cfg,
        }
    }

    /// Returns the resolution of the image
    pub fn get_resolution(&self) -> IVec2 {
        self.get_cfg().resolution
    }

    /// Adds depth of field to camera space ray and transform to world space ray
    fn add_dof(xo_local: Point, wi_local: Direction, cfg: &CameraConfig) -> Ray {
        let (xo_local, wi_local) = if cfg.lens_radius == 0.0 {
            (xo_local, wi_local)
        } else {
            let lens_xy = cfg.lens_radius
                * rand_utils::square_to_disk(rand_utils::unit_square());
            let lens_xyz = lens_xy.extend(0.0);

            let focus_distance = cfg.focal_length / wi_local.z;

            let focus_xyz = focus_distance * wi_local;

            (lens_xyz, focus_xyz - lens_xyz)
        };

        let xo = cfg.camera_to_world.transform_point3(xo_local);
        let wi = cfg.camera_to_world.transform_vector3(wi_local);
        Ray::new(xo, wi)
    }

    /// Generates a ray given a point in raster space `\[0,width\] x \[0,height\]`
    pub fn generate_ray(&self, raster_xy: Vec2) -> Ray {
        match self {
            Self::Perspective(cfg) => {
                let wi_local = cfg.raster_to_camera(raster_xy).normalize();
                Self::add_dof(Point::ZERO, wi_local, cfg)
            }
            Self::Orthographic(cfg, scale) => {
                let xo_local = *scale * cfg.raster_to_camera(raster_xy);
                Self::add_dof(xo_local, Direction::Z, cfg)
            }
        }
    }

    /// Samples a ray leaving from the lens of the camera towards `xi`
    pub fn sample_towards(&self, xi: Point, rand_sq: Vec2) -> Ray {
        let cfg = self.get_cfg();
        let xo_local = rand_utils::square_to_disk(rand_sq).extend(0.0)
            * cfg.lens_radius;
        let xo = cfg.point_to_world(xo_local);

        let wi = (xi - xo).normalize();

        Ray::new(xo, wi)
    }

    /// Probability that `ro` towards `xi` got sampled
    pub fn sample_towards_pdf(&self, ro: &Ray, xi: Point) -> Float {
        let cfg = self.get_cfg();
        let xo = ro.origin;
        let wi = ro.dir;
        let ng = cfg.normal_to_world(Normal::Z);

        let lens_area = if cfg.lens_radius == 0.0 {
            1.0
        } else {
            cfg.lens_radius * cfg.lens_radius * crate::PI
        };

        let pdf = xi.distance_squared(xo) / (ng.dot(wi) * lens_area);
        pdf.max(0.0)
    }

    /// PDF for `wi` direction.
    pub fn pdf_wi(&self, wi: Direction) -> Float {
        let cfg = self.get_cfg();
        let wi_local = cfg.direction_to_local(wi);
        let cos_theta = spherical_utils::cos_theta(wi_local);

        if cos_theta <= 0.0 {
            0.0
        } else {
            let area_coeff = {
                let res = self.get_resolution();
                let res = Vec2::new(
                    res.x as Float,
                    res.y as Float,
                );
                let min_res = res.min_element();
                let screen_bounds = res / min_res;
                screen_bounds.x * screen_bounds.y
            };

            1.0 / (area_coeff * cos_theta * cos_theta * cos_theta)
        }
    }

    /// Incident importance for the ray `ro` starting from the camera lens
    pub fn importance_sample(&self, ro: &Ray) -> FilmSample {
        match self {
            Self::Orthographic(..) => unimplemented!(),
            Self::Perspective(cfg) => {
                let wi = ro.dir;
                let wi_local = cfg.direction_to_local(wi);
                let cos_theta = spherical_utils::cos_theta(wi_local);
                if cos_theta <= 0.0 {
                    return FilmSample::default();
                }

                // compute point in raster space
                let fl = if cfg.lens_radius == 0.0 {
                    1.0 / cos_theta
                } else {
                    cfg.focal_length / cos_theta
                };
                let focus = ro.at(fl);
                let focus_local = cfg.point_to_local(focus);
                let raster_xy = cfg.camera_to_raster(focus_local);

                let lens_area = if cfg.lens_radius == 0.0 {
                    1.0
                } else {
                    crate::PI * cfg.lens_radius * cfg.lens_radius
                };
                let albedo = self.pdf_wi(wi) / (lens_area * cos_theta);

                FilmSample::new(Color::splat(albedo), raster_xy, true)
            }
        }
    }
}
