# Unit Testing Guide for Mass Inflation Code

## Overview

This guide explains the proper numerical regression testing framework for validating changes to the mass inflation simulation code. The key insight is that **we test the physics output** (`data.dat`), not the stdout messages.

## Why Numerical Comparison Instead of `diff`?

Floating-point arithmetic is inherently imprecise due to:
- **Rounding errors**: Different compiler optimizations or CPU architectures may produce slightly different results
- **Order of operations**: `(a + b) + c` may differ from `a + (b + c)` at machine precision
- **Platform differences**: x86 vs ARM, different BLAS libraries, etc.

A simple `diff` will fail even when results are numerically identical to reasonable precision. Instead, we use **tolerance-based comparison**:
```
|test - reference| ≤ atol + rtol × |reference|
```

Where:
- `rtol` = relative tolerance (default: 1e-5 = 0.001%)
- `atol` = absolute tolerance (default: 1e-8)

## Test Framework Components

### 1. Python Comparison Tool (`test/compare_numerical.py`)

A robust numerical comparison tool that:
- Reads multi-column floating-point data
- Skips comment lines (lines starting with `#`)
- Compares with configurable tolerances
- Provides detailed diagnostics when tests fail

**Standalone usage:**
```bash
./test/compare_numerical.py reference_data.dat test_data.dat
./test/compare_numerical.py ref.dat test.dat --rtol 1e-6 --atol 1e-10
./test/compare_numerical.py ref.dat test.dat --columns u v r phi sigma mass
```

### 2. Test Runner Script (`test/run_tests.sh`)

Orchestrates the full test workflow:
1. Builds the code (if needed)
2. Runs the simulation
3. Compares `data.dat` against reference
4. Reports pass/fail with diagnostics

### 3. Makefile Integration

Convenient make targets for testing:
```bash
make test                  # Run regression tests
make test-save-reference   # Save current output as reference
```

## Testing Workflow

### Initial Setup (Before Making Changes)

1. **Ensure code is working correctly** - Review current results, verify physics
2. **Create reference data:**
   ```bash
   make test-save-reference
   ```
   This creates `test/reference_data.dat` from the current code

3. **Commit reference data to version control:**
   ```bash
   git add test/reference_data.dat
   git commit -m "Add reference data for regression tests"
   ```

### Making Changes

1. **Modify the code** (refactor, optimize, modernize, etc.)

2. **Run tests:**
   ```bash
   make test
   ```

3. **Interpret results:**

   **✓ PASS**: Output matches within tolerance
   ```
   ✓ All 1000 rows and 8 columns match within tolerance
   
   ==========================================
     ALL TESTS PASSED
   ==========================================
   ```

   **✗ FAIL**: Numerical differences detected
   ```
   ✗ Numerical differences detected:
   
     Total differences: 250 out of 8000 values
   
     Per-column statistics:
       r         :     50 diffs, max_abs=1.234e-04, max_rel=2.456e-06
       sigma     :    200 diffs, max_abs=5.678e-03, max_rel=1.234e-05
   
     First 10 differences:
        Row    Col       Reference            Test      Abs Diff    Rel Diff
        42      r     5.042008e+00    5.042010e+00     2.000e-06   3.968e-07
   ```

### When Tests Fail

**Non-intentional changes (bugs):**
- Review the code changes
- Use the diagnostic output to identify which physics variables changed
- Debug and fix the issue
- Re-run tests

**Intentional changes (improvements):**
- Verify the new results are physically correct
- Document why results changed
- Update reference data:
  ```bash
  make test-save-reference
  git add test/reference_data.dat
  git commit -m "Update reference: improved integration scheme"
  ```

## Adjusting Tolerances

Default tolerances (rtol=1e-5, atol=1e-8) work for most cases, but you may need to adjust:

**Tighter tolerances** (detecting smaller changes):
```bash
# Edit test/run_tests.sh
RTOL="1e-7"  # More stringent
ATOL="1e-10"
```

**Looser tolerances** (allowing more variation):
```bash
RTOL="1e-4"  # Less stringent
ATOL="1e-6"
```

Choose tolerances based on:
- Physical significance of quantities
- Expected numerical precision of algorithm
- Compiler optimization effects

## Best Practices

### ✓ DO:
- Create reference data from **validated, correct** runs
- Run tests **frequently** during development
- Commit reference data to version control
- Document **why** reference data changes
- Use version control to track test evolution
- Test on **multiple platforms** if portability matters

### ✗ DON'T:
- Blindly update reference without understanding differences
- Use overly loose tolerances that hide bugs
- Test trivial quantities (stdout formatting)
- Ignore systematic errors in specific columns
- Skip tests when "just refactoring"

## Advanced Usage

### Testing Specific Columns

If you only modified code affecting certain physics variables:
```python
# In compare_numerical.py, you can add column filtering
# Or check specific columns manually
```

### Multiple Test Cases

Create different reference files for different physics scenarios:
```bash
test/
  reference_data_schwarzschild.dat
  reference_data_reissner_nordstrom.dat
  reference_data_charged_scalar.dat
```

### Continuous Integration

For automated testing in CI/CD:
```yaml
# .github/workflows/test.yml
- name: Run regression tests
  run: |
    make release
    make test
```

## Column Descriptions

The `data.dat` file contains these physics quantities:

| Column | Variable | Description |
|--------|----------|-------------|
| 1 | u | Retarded null coordinate |
| 2 | v | Advanced null coordinate |
| 3 | r | Radial coordinate |
| 4 | phi | Scalar field |
| 5 | sigma | Conformal factor |
| 6 | mass | Mass function |
| 7 | drdv | Derivative ∂r/∂v |
| 8 | Ricci | Ricci scalar |

Understanding what each column represents helps interpret test failures and determine if differences are physically significant.

## Debugging Failed Tests

1. **Check which columns differ:**
   - Look at "Per-column statistics" in test output
   - Focus on columns with most differences

2. **Examine spatial/temporal patterns:**
   - Are differences only at late times (large v)?
   - Are differences only at large radius?
   - Are differences growing or bounded?

3. **Compare visually:**
   ```bash
   # Plot both datasets
   python3 -c "
   import matplotlib.pyplot as plt
   import numpy as np
   ref = np.loadtxt('test/reference_data.dat')
   test = np.loadtxt('test/output_data.dat')
   plt.plot(ref[:, 0], ref[:, 2], label='reference')
   plt.plot(test[:, 0], test[:, 2], label='test', linestyle='--')
   plt.xlabel('u')
   plt.ylabel('r')
   plt.legend()
   plt.show()
   "
   ```

4. **Use version control to bisect:**
   ```bash
   git bisect start
   git bisect bad              # Current version fails
   git bisect good v1.0        # Known good version
   # Git will checkout commits; run 'make test' at each
   ```

## Debugging with Debug Builds

### Known Issues

**Debug builds may fail with floating-point exceptions**: The debug build enables FPE traps (`-ffpe-trap=invalid,zero,overflow`) to catch numerical errors early. The current code has a latent FPE that occurs during execution. For regression testing, use the `release` build (which doesn't trap FPE), or investigate the FPE using gdb with the debug build.

### Using Debug Builds

The debug build enables runtime checks useful for catching bugs:
```bash
make debug
./bin/mass_inf-debug
```

This enables:
- Array bounds checking (`-fcheck=all`)
- Floating-point exception traps (`-ffpe-trap=invalid,zero,overflow`)
- Runtime warnings (`-Wall -Wextra`)

### Debugging with GDB

```bash
make debug
gdb ./bin/mass_inf-debug
(gdb) run
# When it crashes or traps:
(gdb) backtrace
(gdb) print variable_name
(gdb) info locals
```

Useful GDB commands for Fortran:
- `backtrace` - Show call stack
- `frame N` - Switch to frame N
- `print var` - Print variable value
- `info args` - Show subroutine arguments
- `continue` - Continue execution
- `break function_name` - Set breakpoint

## Migration from Old Tests

The old test system (`make benchmark`) compared stdout output, which is not useful for physics validation. The new system:

| Old System | New System |
|------------|------------|
| Tests stdout formatting | Tests actual physics data |
| Uses simple `diff` | Uses numerical tolerances |
| `test/benchmark.dat` (stdout) | `test/reference_data.dat` (data.dat) |
| `make benchmark` | `make test-save-reference` |
| `make test` (old) | `make test` (new) |

The old targets are still available as `make test-benchmark` for compatibility, but should be considered deprecated.

## Summary

This testing framework provides:
- **Confidence** when refactoring/modernizing code
- **Early detection** of bugs and regressions  
- **Documentation** of how code behavior evolves
- **Foundation** for continuous integration

The key principle: **Test the physics, not the formatting.**
