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

Output lands in a directory named after the parameter file (basename) under the current working directory or under `output_base_dir` if set in `&grid`.

## Namelist layout
- `&grid`: domain bounds, steps (`du`, `dv`), AMR flag and tolerance, output sampling (`output_du`, `output_dv`), progress cadence, optional `output_base_dir`.
- `&model`: model-specific parameters (fields depend on the selected model).

See [examples](examples/) for working templates.

## Models
- RN (4D Reissner-Nordström, default): build with `make` or `make rn`.
- RNLD (D-dimensional Reissner-Nordström-Lambda): build with `make MODEL=rnld`.
- Fixed (scalar field evolution only on fixed RN background): `make MODEL=fixed`.

## Add a model
1) Copy [src/models/template_model.f90](src/models/template_model.f90) to `src/models/<name>_model.f90`.
2) Implement the necessary routines and parameters as specified in [template_model.f90](src/models/template_model.f90) (see also the other implementations in [src/models/](src/models/)).
3) Build with `make MODEL=<name>` or add a convenience target similar to `rn` in the [Makefile](Makefile).

## Tests
- `make tests` runs the regression suite for all models with benchmark data (see the [test](test) folder).
- `make MODEL=<model> test` runs the regression suite for a specific model.
- See [test/README.md](test/README.md) for details.

## Authors

Miguel Zilhão
