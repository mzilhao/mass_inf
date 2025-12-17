# Test Directory

## Quick Start

```bash
# Run all regression cases for default model (rnld)
make test

# Run regression cases for a specific model (looks in test/<model>/)
make MODEL=flat test
```

## Layout

- Per-model folders are required: put configs under `test/<model>/` and references alongside, e.g., `test/flat/D4_config00.nml` with `test/flat/D4_config00/`.
- Each test case has:
  - A config file `*.nml`
  - A reference folder with the same basename, e.g., `D4_config00/`, containing:
    - `fields.dat`       (u, r, phi, sigma)
    - `derivatives.dat`  (u, dr/du, dr/dv)
    - `diagnostics.dat`  (u, mass, Ricci)
- `compare_numerical.py` - Numerical comparator with tolerances
- `run_tests.sh` - Orchestrates builds, runs, and comparisons

## How It Works

1. Builds optimized binary (`make`)
2. For each `*.nml` in this directory:
  - Runs the simulation with working dir `TESTING/` (outputs land in `TESTING/<case>/`)
  - Compares `fields.dat`, `derivatives.dat`, `diagnostics.dat` against the reference set
  - Uses tolerances: `rtol=1e-5`, `atol=2e-8`
3. Reports pass/fail with per-case details

## Legacy Files (Deprecated)

- `benchmark.dat` - Old stdout-based tests (use `make test-benchmark`)
- `output.dat` - Old test output

The new system tests actual physics output, not stdout formatting.

See `../TESTING.md` for comprehensive documentation.
