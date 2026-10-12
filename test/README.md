# Test Directory


## Layout

- Per-model folders are required under `test/<model>/`
- Each test case has:
  - A config file `*.nml`
  - A reference folder with the same basename, e.g., `RN_test00/`, containing the benchmark data files
- `compare_numerical.py` - Numerical comparator with tolerances
- `run_tests.sh` - Orchestrates builds, runs, and comparisons

## How It Works

1. The model binary is built first by the Makefile (`make [MODEL=<model>] test`)
2. For each `*.nml` in the testing directory:
  - Runs the simulation with working dir `TESTING/<model>/` (outputs land in `TESTING/<model>/<case>/`)
  - Compares output `*.dat` files against the reference set with the set tolerances
3. Reports pass/fail with per-case details

## Notes
- The model directory name under `test/` must match the `MODEL` spelling exactly, including case (e.g. `test/Bardeen` for `make MODEL=Bardeen`).
- `make tests` aborts at the first failing model.
- To regenerate reference data after an intentional output change: run each affected case from a scratch directory (the binary writes into a folder named after the parfile basename, under the CWD or `output_base_dir`) and copy the new `*.dat` files into `test/<model>/<case>/`.
