FC = gfortran
FCFLAGS_RELEASE = -O3 -march=native -fno-unsafe-math-optimizations
FCFLAGS_DEBUG = -g -Wall -Wextra -Wno-unused-dummy-argument -fcheck=all \
                -ffpe-trap=invalid,zero,overflow -Wno-maybe-uninitialized
MODFLAGS_RELEASE = -J$(OBJDIR_RELEASE)
MODFLAGS_DEBUG = -J$(OBJDIR_DEBUG)
LDFLAGS =

# Directories
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
.PHONY: all release debug clean distclean test test-save-reference help

# Default target: release build
all: release

# Release build (optimized)
release: $(BINDIR)/mass_inf

# Debug build (with bounds checking and traps)
debug: $(BINDIR)/mass_inf-debug

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

# Clean build artifacts
clean:
	rm -rf $(OBJDIR) $(BINDIR) *.mod *.MOD

# Full cleanup (removes all generated files)
distclean: clean
	rm -f data.dat test/output_data.dat

# Run numerical regression tests
test: release
	@./test/run_tests.sh

# Create/update reference data for regression tests
test-save-reference: release
	@echo "Saving current output as reference data..."
	@$(BINDIR)/mass_inf > /dev/null 2>&1
	@cp data.dat test/reference_data.dat
	@echo "✓ Reference data saved to test/reference_data.dat"

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
	@echo "  make help                 - Show this message"

