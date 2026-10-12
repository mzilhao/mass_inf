# AGENTS.md

Fortran (gfortran) double-null mass-inflation solver in spherical symmetry.
Each executable is built with **one** physics model compiled in; a run reads a
single namelist file containing a `&grid` block (numerics) and a `&model`
block (physics, defined per model).

## Build
- `make` builds the default `rn` model -> `bin/mass_inf-rn`. Others: `make MODEL=rnld|fixed|Bardeen` (also `make rn`); `make models` builds rn+rnld+fixed.
- Model is a compile-time choice: `src/models/<name>_model.f90` is compiled as `model_mod.o` and linked in. Binaries are always model-suffixed — the Makefile never produces an unsuffixed `mass_inf`. Switching `MODEL` always recompiles the model object (forced via a `FORCE` pseudo-target), so a stale binary from a previously built model cannot be silently relinked.
- `Bardeen` has no convenience target (build with `make MODEL=Bardeen`), but it is covered by `make tests` (bench data under `test/Bardeen/`). `rn` selects initial data via `ics` in `&model`: `'gaussian'` (default, compact pulse) or `'log'` (the former noncomp model; old noncomp parfiles need `ics = 'log'` added).
- `make [MODEL=<m>] debug` -> `bin/mass_inf-<m>-debug` with `-fcheck=all -ffpe-trap=invalid,zero,overflow` (use when chasing NaNs/crashes).
- Release `FCFLAGS` uses `-march=native` — binaries are machine-specific.

### Module dependencies are manual
`src/*.f90` files are collected by wildcard, but `.o` module dependencies are
listed by hand under "Module dependencies" in the Makefile. When adding a
top-level module, add its dependency lines there too.

## Run
- `./bin/mass_inf-<model> [-f] <config.nml>`; templates in `examples/`.
- Output dir = `output_dir` from `&grid` (default: parfile basename without
  extension), placed under `output_base_dir` (default: CWD; examples use
  `./runs/`).
- The run **fails** if the output directory already exists unless `-f` is given
  (`-f` removes it and copies the parfile inside).

## Tests
- `make tests` — regression suite for every model with data under `test/<model>/`.
  `make [MODEL=<m>] test` for a single model (builds the binary first).
- Layout: `test/<model>/<case>.nml` plus a sibling reference dir `<case>/` of
  `.dat` files. `<model>` must match the MODEL spelling exactly, including
  case (`test/Bardeen` for the Bardeen model). New dirs under `test/` are
  picked up automatically by `make tests`, which now aborts on the first
  failing model. Runs execute in `TESTING/<model>/<case>/` (gitignored).
- Comparisons use `test/compare_numerical.py` with tolerances hardcoded in
  `test/run_tests.sh` (RTOL=1e-5, ATOL=2e-8). Requires python3 + numpy (the
  script activates a conda env `numpy` if present).
- Intentional output changes: regenerate/update the reference `.dat` files in
  the model's test dir.
- Output `.dat` files must keep their `# Columns: ...` header comments — both
  the comparator and the analysis code discover columns from them.

## Adding a model
1. Copy `src/models/template_model.f90` to `src/models/<name>_model.f90`
   (filename matters: it is compiled as `model_mod.o`). Optionally add a
   convenience target like `rn` in the Makefile.
2. Implement the stable API: `model_config_mod` (type `model_config` + `load`
   reading `&model`) and `model_mod` (`NEQ`, `F`, `init_cond`,
   `open_output_files`, `write_output`, `write_constraints`,
   `close_output_files`).
3. To be covered by `make tests`, add configs + reference dirs under
   `test/<model>/`.

## Code layout
- One module per file in `src/`: `precision` (kind `dp`), `grid_config`
  (parses `&grid`, derives `Nu`/`Nv`/`Nu_max`), `pde_stepper`/`evolve` (RHS
  integrator, generic in `NEQ`, Picard iteration), `AMR`, `utils`,
  `polint`; `mass_inf.f90` is the main program. Model code lives only in
  `src/models/<name>_model.f90`.
- All models store fields as `(r, sigma, phi)` in `h()`. Keep new models in
  this order if possible.
- `rnld`'s `lambda` is not the bare Einstein Λ for D≠4: the metric function is
  `f = 1 − 2 m0/r^(D−3) + q²/r^(2(D−3)) − λ r²/3`, i.e.
  `λ = 6·Λ_Einstein/((D−1)(D−2))` (they coincide only at D=4).
- AMR refines the u-direction grid at runtime using `plus`/`minus` linked
  lists; arrays are sized `Nu_max = 2*(u_max-u_min)/reldiff_max` — beware
  memory with small `reldiff_max`.

## Analysis
- `analysis/read_data.py` (importable package `analysis`): `MassInflationData(folder)`
  reads ASCII run outputs; fields discovered from `# Columns:` headers of the
  `.dat` files. Simulation output goes to `runs/` (gitignored).
- `analysis/test_analysis.py` is a `# %%` cell-style scratch/plot script
  (needs matplotlib).