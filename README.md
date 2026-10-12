# mass_inf

Double null mass inflation solver (in spherical symmetry) with pluggable physics model implementations. The executable is built with one model at a time and reads a single namelist file with two blocks: `&grid` for numerics and `&model` for physics/model parameters.


## Requirements
- gfortran
- POSIX shell tools (for `make`, tests)
- python, numpy (for tests and analysis scripts)

## Build
- Default (RN model):
   - `make` → `bin/mass_inf-rn`
- Other models:
   - `make MODEL=<model>` → `bin/mass_inf-<model>`

## Run
Use a namelist file containing both `&grid` and `&model` blocks. Example:

```zsh
./bin/mass_inf-rn examples/RN_config00.nml
```

Output lands in a directory named after the parameter file (basename) under the current working directory or under `output_base_dir` if set in `&grid`. The run aborts if the output directory already exists; pass `-f` to remove it and overwrite (e.g. `./bin/mass_inf-rn -f examples/RN_config00.nml`).

## Namelist layout
- `&grid`: domain bounds, steps (`du`, `dv`), AMR flag and tolerance, output sampling (`output_du`, `output_dv`), progress cadence, optional `output_base_dir`.
- `&model`: model-specific parameters (fields depend on the selected model).

See [examples](examples/) for working templates.

## Models
- RN (4D Reissner-Nordström, default): build with `make` or `make rn`. The scalar-field initial data on the `u = u_min` slice is selected with `ics` in `&model`: `'gaussian'` (default; a compact pulse) or `'log'` (`phi = A log(v/v_min)`, a permanently incoming flux — the former `noncomp` model; old parfiles need `ics = 'log'` added).
- Bardeen (regular Bardeen black hole sourced by nonlinear electrodynamics): `make MODEL=Bardeen`.
- RNLD (D-dimensional charged black hole with cosmological constant): `make MODEL=rnld`. Conventions: the metric function contains `q` as `q²/r^(2(D−3))` and the cosmological term as `− lambda r²/3`; `lambda` coincides with the bare Einstein `Λ` only for `D = 4` (generally `lambda = 6Λ/((D−1)(D−2))`).
- Fixed (the geometry is prescribed and never evolved: `r` and `sigma` remain on their initial slices, only the scalar field evolves): `make MODEL=fixed`.

With `A = 0`, the RN, Bardeen and RNLD models reproduce their exact background solutions — a good self-check is that the `mass` diagnostic stays constant. `make models` builds the rn, rnld and fixed binaries.

## Add a model
1) Copy [src/models/template_model.f90](src/models/template_model.f90) to `src/models/<name>_model.f90`.
2) Implement the necessary routines and parameters as specified in [template_model.f90](src/models/template_model.f90) (see also the other implementations in [src/models/](src/models/)).
3) Build with `make MODEL=<name>` or add a convenience target similar to `rn` in the [Makefile](Makefile).

Notes:
- `.o` module dependencies are listed by hand in the "Module dependencies" section of the [Makefile](Makefile); add lines for any new top-level module in `src/`.
- All models store the field arrays in the order `(r, sigma, phi)`.
- Regression cases live under `test/<model>/` (see [Tests](#tests)).

## Tests
- `make tests` runs the regression suite for all models with benchmark data (see the [test](test) folder). Model directories must match the `MODEL` spelling exactly, including case (e.g. `test/Bardeen`), and the suite aborts at the first failing model.
- `make MODEL=<model> test` runs the regression suite for a specific model.
- If a change intentionally alters the output, regenerate the reference data: run each affected case from any scratch directory and copy the produced `*.dat` files into `test/<model>/<case>/`.
- Tolerances (RTOL = 1e-5, ATOL = 2e-8) are defined in [test/run_tests.sh](test/run_tests.sh) and can be overridden from the environment.
- See [test/README.md](test/README.md) for details.

## Authors

Miguel Zilhão  
Ana Bokulić
