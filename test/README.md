# Test Directory

## Quick Start

```bash
# Run tests (compares data.dat against reference)
make test

# Save current output as new reference (after validating changes)
make test-save-reference
```

## Files

- `reference_data.dat` - Golden reference data (120,000 rows × 8 columns)
- `compare_numerical.py` - Python tool for numerical comparison with tolerances
- `compare_numerical.awk` - AWK fallback (if numpy unavailable)
- `run_tests.sh` - Test orchestration script

## How It Works

1. Builds optimized binary (`make release`)
2. Runs simulation → produces `data.dat`
3. Compares against `reference_data.dat` using:
   - Relative tolerance: 1e-5 (0.001%)
   - Absolute tolerance: 1e-8
4. Reports pass/fail with detailed diagnostics

## Legacy Files (Deprecated)

- `benchmark.dat` - Old stdout-based tests (use `make test-benchmark`)
- `output.dat` - Old test output

The new system tests actual physics output, not stdout formatting.

See `../TESTING.md` for comprehensive documentation.
