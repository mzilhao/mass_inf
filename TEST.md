# Regression Testing Guide

## Overview
This project includes a simple regression testing framework to validate changes against a saved benchmark output.

## Known Issues

**Debug builds may fail with floating-point exceptions**: The debug build enables FPE traps (`-ffpe-trap=invalid,zero,overflow`) to catch numerical errors early. The current code has a latent FPE that occurs during execution. This is intentional—it helps identify issues for future fixes. For regression testing, use the `release` build (which doesn't trap FPE), or investigate the FPE using gdb with the debug build.

## Workflow

### Initial Setup (Baseline)
Before making changes, establish a golden reference output:

```bash
make benchmark
```

This runs the compiled binary and saves output to `test/benchmark.dat`. Commit this file to version control.

### Running Tests
After making code changes, verify correctness by comparing against the benchmark:

```bash
make test
```

This will:
1. Rebuild with release flags (`make release`) — optimized but without FPE traps
2. Run the binary and capture output to `test/output.dat`
3. Compare `test/output.dat` against `test/benchmark.dat`
4. Report pass/fail

**Pass example:**
```
✓ Test passed: output matches benchmark
```

**Fail example:**
```
✗ Test failed: output differs from benchmark
Differences:
[diff output showing first 20 lines of divergence]
```

### Accepting Changes
If changes are intentional and validated (e.g., algorithm improvement), update the benchmark:

```bash
make benchmark
```

Then commit `test/benchmark.dat`.

## Build Targets

| Target       | Purpose |
|--------------|---------|
| `make all` or `make release` | Build optimized version (default) |
| `make debug` | Build with bounds checking & FPE traps |
| `make test` | Run regression tests (rebuilds debug) |
| `make benchmark` | Generate/update golden output |
| `make clean` | Remove build artifacts (preserve benchmarks) |
| `make distclean` | Full cleanup including benchmarks |
| `make help` | Show all targets |

## Debugging Failed Tests

If tests fail unexpectedly:

1. **View differences:**
   ```bash
   diff test/benchmark.dat test/output.dat | head -50
   ```

2. **Run with debug build:**
   ```bash
   make debug
   ./bin/mass_inf > /tmp/debug.dat 2>&1
   ```
   This enables bounds checking, overflow traps, and warnings.

3. **Compare intermediate values:**
   Modify output in `mass_inf.f90` to log intermediate state if needed.

4. **Check for numerical drift:**
   Small differences may be acceptable if within floating-point precision. Review the diff carefully.

## Notes

- The benchmark is stored in `test/benchmark.dat` and should be version-controlled.
- Output files use double precision, so small numerical variations (ULP level) are normal.
- For large simulations, consider adding a tolerance-based comparison rather than exact matching.
