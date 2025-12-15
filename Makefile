FC = gfortran
FCFLAGS = -O3 -march=native -fno-unsafe-math-optimizations
FCFLAGS_DEBUG = -g -Wall -Wextra -Wno-unused-dummy-argument -fcheck=all \
                -ffpe-trap=invalid,zero,overflow -Wno-maybe-uninitialized
FCFLAGS_PROFILE = -O2 -pg
LDFLAGS =

# Directories
SRCDIR = src
BINDIR = bin
OBJDIR = obj

# Model selection (compile-time)
MODEL ?= rnld
MODELDIR = $(SRCDIR)/models
MODEL_SRC = $(MODELDIR)/$(MODEL)_model_impl.f90
MODEL_OBJ = $(OBJDIR)/model_mod.o
MODEL_OBJ_DEBUG = $(OBJDIR)/debug_model_mod.o

# Binary naming
BIN_BASENAME ?= mass_inf
BIN_SUFFIX ?=
BINARY = $(BINDIR)/$(BIN_BASENAME)$(BIN_SUFFIX)

# Source files and objects (top-level src only)
SOURCES = $(wildcard $(SRCDIR)/*.f90)
OBJECTS = $(patsubst $(SRCDIR)/%.f90,$(OBJDIR)/%.o,$(SOURCES)) $(MODEL_OBJ)
OBJECTS_DEBUG = $(patsubst $(SRCDIR)/%.f90,$(OBJDIR)/debug_%.o,$(SOURCES)) $(MODEL_OBJ_DEBUG)

# Phony targets
.PHONY: all debug profile clean distclean test test-save-reference help rnld dummy template models

# Default target
all: $(BINARY)

# Main executable
$(BINARY): $(OBJECTS) | $(BINDIR)
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

# Debug build
debug: $(BINDIR)/mass_inf-debug

$(BINDIR)/mass_inf-debug: $(OBJECTS_DEBUG) | $(BINDIR)
	$(FC) $(FCFLAGS_DEBUG) -o $@ $^ $(LDFLAGS)

# Profile build
profile: clean
	$(MAKE) FCFLAGS="$(FCFLAGS_PROFILE)" LDFLAGS="-pg" $(BINDIR)/mass_inf-prof

$(BINDIR)/mass_inf-prof: $(OBJECTS) | $(BINDIR)
	$(FC) $(FCFLAGS_PROFILE) -o $@ $^ -pg

# Compile release objects (top-level sources)
$(OBJDIR)/%.o: $(SRCDIR)/%.f90 | $(OBJDIR)
	$(FC) $(FCFLAGS) -J$(OBJDIR) -c -o $@ $<

# Compile debug objects (top-level sources)
$(OBJDIR)/debug_%.o: $(SRCDIR)/%.f90 | $(OBJDIR)
	$(FC) $(FCFLAGS_DEBUG) -J$(OBJDIR) -c -o $@ $<

# Compile selected model implementation
$(MODEL_OBJ): $(MODEL_SRC) | $(OBJDIR)
	$(FC) $(FCFLAGS) -J$(OBJDIR) -c -o $@ $<

$(MODEL_OBJ_DEBUG): $(MODEL_SRC) | $(OBJDIR)
	$(FC) $(FCFLAGS_DEBUG) -J$(OBJDIR) -c -o $@ $<

# Module dependencies
$(OBJDIR)/grid_config.o: $(OBJDIR)/precision.o
$(MODEL_OBJ): $(OBJDIR)/precision.o $(OBJDIR)/grid_config.o
$(OBJDIR)/pde_stepper.o: $(OBJDIR)/precision.o $(MODEL_OBJ)
$(OBJDIR)/evolve.o: $(OBJDIR)/precision.o $(OBJDIR)/pde_stepper.o $(MODEL_OBJ)
$(OBJDIR)/mass_inf.o: $(OBJDIR)/precision.o $(MODEL_OBJ) $(OBJDIR)/evolve.o $(OBJDIR)/grid_config.o $(OBJDIR)/polint.o $(OBJDIR)/utils.o $(OBJDIR)/AMR.o
$(OBJDIR)/utils.o: $(OBJDIR)/precision.o
$(OBJDIR)/polint.o: $(OBJDIR)/precision.o
$(OBJDIR)/AMR.o: $(OBJDIR)/precision.o $(OBJDIR)/grid_config.o $(OBJDIR)/polint.o $(OBJDIR)/utils.o $(MODEL_OBJ)

$(MODEL_OBJ_DEBUG): $(OBJDIR)/debug_grid_config.o
$(OBJDIR)/debug_pde_stepper.o: $(MODEL_OBJ_DEBUG)
$(OBJDIR)/debug_evolve.o: $(OBJDIR)/debug_pde_stepper.o $(MODEL_OBJ_DEBUG)
$(OBJDIR)/debug_mass_inf.o: $(MODEL_OBJ_DEBUG) $(OBJDIR)/debug_evolve.o $(OBJDIR)/debug_grid_config.o $(OBJDIR)/debug_polint.o $(OBJDIR)/debug_utils.o $(OBJDIR)/debug_AMR.o

# Create directories
$(BINDIR) $(OBJDIR):
	mkdir -p $@

# Clean build artifacts
clean:
	rm -rf $(OBJDIR) $(BINDIR) *.mod *.MOD gmon.out

# Full cleanup
distclean: clean
	rm -f data.dat timing.log gprof.txt
	rm -rf TESTING

# Run tests
test: all
	@./test/run_tests.sh

# Convenience targets to build specific models with distinct binaries
rnld:
	$(MAKE) MODEL=rnld BIN_SUFFIX=-rnld all

dummy:
	$(MAKE) MODEL=dummy BIN_SUFFIX=-dummy all

template:
	$(MAKE) MODEL=template BIN_SUFFIX=-template all

# Build all provided models (cleans between builds to avoid cross-contamination)
models:
	$(MAKE) clean
	$(MAKE) rnld
	$(MAKE) clean
	$(MAKE) dummy
	$(MAKE) clean
	$(MAKE) template

# Save reference data
test-save-reference: all
	@echo "Saving current outputs as references for all test cases..."
	@bash test/run_tests.sh --save-reference

# Help
help:
	@echo "Available targets:"
	@echo "  make           - Build optimized binary (default)"
	@echo "  make debug     - Build with debug flags and checks"
	@echo "  make profile   - Build with gprof instrumentation"
	@echo "  make test      - Run regression tests"
	@echo "  make rnld      - Build RNLD model -> bin/mass_inf-rnld"
	@echo "  make dummy     - Build dummy model -> bin/mass_inf-dummy"
	@echo "  make template  - Build template model -> bin/mass_inf-template"
	@echo "  make models    - Build all models (rnld, dummy, template)"
	@echo "  make clean     - Remove build artifacts"
	@echo "  make distclean - Remove all generated files"
	@echo "  make help      - Show this message"
