# Quick Reference: Build & Test

## Standard Workflow

```bash
# Initial setup: create baseline
make benchmark

# Make your changes...
# Edit functions.f90, evolve.f90, etc.

# Run tests to verify
make test

# If output differs, investigate
diff test/benchmark.dat test/output.dat | head -50

# If changes are correct, update benchmark
make benchmark
```

## Build Targets

```bash
make release      # Optimized build (default for 'make all')
make debug        # Build with checks & FPE traps for debugging
make test         # Regression test (compares against benchmark)
make benchmark    # Capture current run as golden reference
make clean        # Remove objects, modules, binaries (keep benchmarks)
make distclean    # Full cleanup including benchmarks
make help         # Show this help text
```

## Debug a Floating-Point Exception

```bash
# 1. Build with debug flags
make debug

# 2. Run under debugger
gdb --args ./bin/mass_inf

# 3. Inside gdb:
(gdb) run
# ...will stop at FPE with backtrace

# 4. Examine variables around the crash
(gdb) print variable_name
(gdb) bt  # backtrace
```

## Verify Code After Edits

```bash
# Quick check
make clean && make release && make test

# Comprehensive check (bounds checking, FPE traps)
make debug
./bin/mass_inf > /tmp/test.dat 2>&1

# Compare outputs
diff test/benchmark.dat /tmp/test.dat | less
```

## File Organization

```
project/
├── Makefile              # Build system
├── *.f90                 # Source files
├── bin/                  # Compiled binaries (git-ignored)
├── obj/                  # Object files & modules (git-ignored)
├── test/
│   ├── benchmark.dat     # Golden reference (COMMITTED)
│   └── output.dat        # Latest test run (git-ignored)
├── TEST.md               # Testing guide
└── IMPROVEMENTS.md       # What was improved
```

## Common Commands

| Command | Purpose |
|---------|---------|
| `make test` | Most common: test your changes |
| `make benchmark` | Save current as baseline after validation |
| `make debug && ./bin/mass_inf` | Debug with FPE detection |
| `make clean && make release` | Full rebuild from scratch |
| `make help` | Show all targets |
