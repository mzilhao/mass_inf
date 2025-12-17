# mass_inf

Mass inflation solver with pluggable model implementations. The executalble is built with one model at a time and reads a single namelist file with two blocks: `&grid` for numerics and `&model` for physics/model parameters.


## Requirements
- gfortran
- POSIX shell tools (for `make`, tests)

## Build
- Default (RNLD model):
   - `make` → `bin/mass_inf-rnld`
- Other models:
   - `make flat` → `bin/mass_inf-flat`
   - `make template` (scaffolding only)
- Debug/profile:
   - `make rnld-debug`, `make flat-debug`
   - `make rnld-profile`, `make flat-profile`
- Clean: `make clean`

## Run
Use a namelist file containing both `&grid` and `&model` blocks. Example:

```zsh
./bin/mass_inf-rnld examples/D4_config00.nml
```

Output lands in a directory named after the parameter file (basename) under the current working directory or under `output_base_dir` if set in `&grid`.

## Namelist layout
- `&grid`: domain bounds, steps (`du`, `dv`), AMR flag and tolerance, output sampling (`output_du`, `output_dv`), progress cadence, optional `output_base_dir`.
- `&model`: model-specific parameters (fields depend on the selected model).

See [examples/D4_config00.nml](examples/D4_config00.nml) for a working template.

## Models
- RNLD (default): build with `make` or `make rnld`.
- Flat: `make flat`.
- Template scaffold: `make template` (for reference only).

## Add a model
1) Copy [src/models/template_model.f90](src/models/template_model.f90) to `src/models/<name>_model.f90`.
2) Implement `model_config_mod` and `model_mod` (set `NEQ`, provide `load`, `F`, `init_cond`, output hooks).
3) Build with `make MODEL=<name> BIN_SUFFIX=-<name>` or add a convenience target similar to `rnld`/`flat` in the [Makefile](Makefile).

## Tests
- `make test` runs the regression suite (uses the default model build) and compares outputs in [test](test).


## Authors

Miguel Zilhão

