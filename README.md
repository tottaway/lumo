## Lumo
[![crates.io](https://img.shields.io/crates/v/lumo)](https://crates.io/crates/lumo)
[![docs.rs](https://img.shields.io/docsrs/lumo)](https://docs.rs/lumo)
[![Coverage](https://img.shields.io/coverallsCoverage/github/ekarpp/lumo)](https://coveralls.io/github/ekarpp/lumo)

Lumo is a CPU based multithreaded rendering engine. Made with the goal of learning Rust and physically based rendering :)

### Features
* Path tracing and bidirectional path tracing with [MIS](http://iliyan.com/publications/ImplementingVCM)
* [Cook-Torrance microfacet BSDF](https://doi.org/10.1145/357290.357293) with [Beckmann and GGX](http://dx.doi.org/10.2312/EGWR/EGSR07/195-206)
* .obj and .mtl file parsing
* [Surface area hierarchy based kD-trees](https://www.irisa.fr/prive/kadi/Sujets_CTR/kadi/Kadi_sujet2_article_Kdtree.pdf)

### Gallery
![Bust of Nefertiti](https://i.imgur.com/XuLT7Wy.png)
![Cornell box](https://i.imgur.com/0EozvDq.png)
![Conference room](https://i.imgur.com/7YyWNKr.png)
![Circle of spheres](https://i.imgur.com/zraIbaH.png)

### Usage
Once the repository is cloned, the `examples/` folder contains scenes, that can be ran with `cargo`:

```bash
cargo run --example hello_sphere
```

The renderer can be configured either through its setter methods in the examples or partially through the CLI:

```
Usage: hello_sphere [-s <samples>] [-t <threads>] [-d] [-b]

Optional CLI configuration of renderer. Renderer setter methods have priority.

Options:
  -s, --samples     number of samples per pixel (defaults to 1)
  -t, --threads     number of threads used (defaults to all)
  -d, --direct      use direct light integrator instead of path tracing
  -b, --bdpt        use bidirectional path tracing instead of path tracing
  --help            display usage information
```

### References
* [Physically Based Rendering](https://www.pbr-book.org/)
* [Moving Frostbite to Physically Based Rendering](https://seblagarde.files.wordpress.com/2015/07/course_notes_moving_frostbite_to_pbr_v32.pdf)
* [Eric Veach's PhD Thesis](http://graphics.stanford.edu/papers/veach_thesis/)
