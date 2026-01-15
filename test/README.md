# Test Directory


## Layout

- Per-model folders are required under `test/<model>/`
- Each test case has:
  - A config file `*.nml`
  - A reference folder with the same basename, e.g., `RN_test00/`, containing the benchmark data files
- `compare_numerical.py` - Numerical comparator with tolerances
- `run_tests.sh` - Orchestrates builds, runs, and comparisons

## How It Works

1. Builds binary (`make MODEL=<model>`) if not already built
2. For each `*.nml` in the testing directory:
  - Runs the simulation with working dir `TESTING/` (outputs land in `TESTING/<case>/`)
  - Compares output `*.dat` files against the reference set with the set tolerances
3. Reports pass/fail with per-case details
