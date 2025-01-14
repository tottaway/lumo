use crate::tracer::Integrator;
use clap::Parser;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
/// Optional CLI configuration of renderer. Renderer setter methods have priority.
pub struct TracerCli {
    /// number of samples per pixel (defaults to 1)
    #[arg(short, long, default_value_t = 1)]
    pub samples: i32,

    /// number of threads used (defaults to all)
    #[arg(short, long)]
    pub threads: Option<usize>,

    /// use direct light integrator instead of path tracing
    #[arg(short, long = "direct")]
    pub direct_light: bool,

    /// use bidirectional path tracing instead of path tracing
    #[arg(short, long = "bdpt")]
    pub bd_path_trace: bool,
}

impl TracerCli {
    /// Sets the configured number of threads. Called by renderer on creation.
    pub fn set_threads(&self) {
        if let Some(t) = self.threads {
            rayon::ThreadPoolBuilder::new()
                .num_threads(t)
                .build_global()
                .unwrap();
        }
    }

    /// Get the configured integrator.
    pub fn get_integrator(&self) -> Integrator {
        if self.direct_light {
            Integrator::DirectLight
        } else if self.bd_path_trace {
            Integrator::BDPathTrace
        } else {
            Integrator::PathTrace
        }
    }
}
