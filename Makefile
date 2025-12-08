FC = gfortran
FCFLAGS_RELEASE = -O3 -march=native -fno-unsafe-math-optimizations
FCFLAGS_DEBUG = -g -Wall -Wextra -Wno-unused-dummy-argument -fcheck=all \
                -ffpe-trap=invalid,zero,overflow -Wno-maybe-uninitialized
# Compiler flags to put .mod files in the object directory
MODFLAGS_RELEASE = -J$(OBJDIR_RELEASE)
MODFLAGS_DEBUG = -J$(OBJDIR_DEBUG)
LDFLAGS =

# List of executables to be built
PROGRAMS = mass_inf
SRCDIR = .
BINDIR = bin
OBJDIR = obj
OBJDIR_RELEASE = $(OBJDIR)/release
OBJDIR_DEBUG = $(OBJDIR)/debug

# Create directories if they don't exist
$(shell mkdir -p $(BINDIR) $(OBJDIR_RELEASE) $(OBJDIR_DEBUG))

# Automatically discover source files and compute object files
SOURCES = $(wildcard $(SRCDIR)/*.f90)
OBJECTS_RELEASE = $(patsubst $(SRCDIR)/%.f90,$(OBJDIR_RELEASE)/%.o,$(SOURCES))
OBJECTS_DEBUG = $(patsubst $(SRCDIR)/%.f90,$(OBJDIR_DEBUG)/%.o,$(SOURCES))

# Phony targets
.PHONY: all debug release clean test help

# Default target: release build
all: release

# Release build (optimized)
release: FCFLAGS = $(FCFLAGS_RELEASE)
release: $(addprefix $(BINDIR)/,$(PROGRAMS))

# Debug build (with bounds checking and traps)
debug: FCFLAGS = $(FCFLAGS_DEBUG)
debug: $(addprefix $(BINDIR)/,$(PROGRAMS))-debug

# Link release executable
$(BINDIR)/mass_inf: $(OBJECTS_RELEASE)
	$(FC) $(FCFLAGS_RELEASE) -o $@ $^ $(LDFLAGS)

# Link debug executable
$(BINDIR)/mass_inf-debug: $(OBJECTS_DEBUG)
	$(FC) $(FCFLAGS_DEBUG) -o $@ $^ $(LDFLAGS)

# Compile release objects
$(OBJDIR_RELEASE)/%.o: $(SRCDIR)/%.f90
	$(FC) $(FCFLAGS_RELEASE) $(MODFLAGS_RELEASE) -c -o $@ $<

# Compile debug objects
$(OBJDIR_DEBUG)/%.o: $(SRCDIR)/%.f90
	$(FC) $(FCFLAGS_DEBUG) $(MODFLAGS_DEBUG) -c -o $@ $<

# Module dependencies for release build
# functions.f90 contains physics_config_mod, must be built first
$(OBJDIR_RELEASE)/pde_stepper.o: $(OBJDIR_RELEASE)/functions.o
$(OBJDIR_RELEASE)/evolve.o: $(OBJDIR_RELEASE)/pde_stepper.o $(OBJDIR_RELEASE)/functions.o
$(OBJDIR_RELEASE)/mass_inf.o: $(OBJDIR_RELEASE)/functions.o $(OBJDIR_RELEASE)/evolve.o

# Module dependencies for debug build
$(OBJDIR_DEBUG)/pde_stepper.o: $(OBJDIR_DEBUG)/functions.o
$(OBJDIR_DEBUG)/evolve.o: $(OBJDIR_DEBUG)/pde_stepper.o $(OBJDIR_DEBUG)/functions.o
$(OBJDIR_DEBUG)/mass_inf.o: $(OBJDIR_DEBUG)/functions.o $(OBJDIR_DEBUG)/evolve.o

# Clean build artifacts (keep benchmark data)
clean:
	rm -rf $(OBJDIR) $(BINDIR) *.mod *.MOD

# Full cleanup (removes all generated files including benchmarks)
distclean: clean
	rm -f test/*.out test/*.dat

# Run numerical regression tests (new proper testing framework)
test: release
	@echo "Running numerical regression tests..."
	@./test/run_tests.sh

# Create/update reference data for regression tests
test-save-reference: release
	@echo "Saving current output as reference data..."
	@mkdir -p test
	@$(BINDIR)/mass_inf > /dev/null 2>&1
	@cp data.dat test/reference_data.dat
	@echo "✓ Reference data saved to test/reference_data.dat"
	@echo "  This will be used for future regression testing."

# Run old stdout-based benchmark tests (deprecated, kept for compatibility)
test-benchmark: release
	@echo "Running old benchmark tests (deprecated)..."
	@mkdir -p test
	@$(BINDIR)/mass_inf > test/output.dat 2>&1
	@if [ -f test/benchmark.dat ]; then \
		if diff -q test/output.dat test/benchmark.dat > /dev/null 2>&1; then \
			echo "✓ Test passed: output matches benchmark"; \
		else \
			echo "✗ Test failed: output differs from benchmark"; \
			echo "Differences:"; \
			diff test/output.dat test/benchmark.dat | head -20; \
		fi; \
	else \
		echo "⚠ No benchmark found. Saving current output as benchmark..."; \
		cp test/output.dat test/benchmark.dat; \
	fi

# Save current stdout as benchmark (deprecated)
benchmark: release
	@echo "Creating benchmark from current run..."
	@mkdir -p test
	@$(BINDIR)/mass_inf > test/benchmark.dat 2>&1
	@echo "Benchmark saved to test/benchmark.dat"

# Help target
help:
	@echo "Available targets:"
	@echo "  make all                  - Build release (optimized) version (default)"
	@echo "  make release              - Build optimized version"
	@echo "  make debug                - Build debug version with checks"
	@echo "  make test                 - Run numerical regression tests"
	@echo "  make test-save-reference  - Save current data.dat as reference"
	@echo "  make clean                - Remove build artifacts"
	@echo "  make distclean            - Remove all generated files"
	@echo "  make release      - Explicitly build release version (bin/mass_inf)"
	@echo "  make debug        - Build debug version with checks (bin/mass_inf-debug)"
	@echo "  make test         - Run regression tests against benchmark"
	@echo "  make benchmark    - Save current output as benchmark"
	@echo "  make clean        - Remove build artifacts (keep benchmarks)"
	@echo "  make distclean    - Remove all generated files including benchmarks"
	@echo "  make help         - Show this message"

