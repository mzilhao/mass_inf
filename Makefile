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
MODEL_SRC = $(MODELDIR)/$(MODEL)_model.f90
MODEL_OBJ = $(OBJDIR)/model_mod.o
MODEL_OBJ_DEBUG = $(OBJDIR)/debug_model_mod.o

# Binary naming
BIN_BASENAME ?= mass_inf
# Always suffix binaries with the selected model to avoid unsuffixed 'mass_inf'
BIN_SUFFIX ?= -$(MODEL)
BINARY = $(BINDIR)/$(BIN_BASENAME)$(BIN_SUFFIX)

# Source files and objects (top-level src only)
SOURCES = $(wildcard $(SRCDIR)/*.f90)
OBJECTS = $(patsubst $(SRCDIR)/%.f90,$(OBJDIR)/%.o,$(SOURCES)) $(MODEL_OBJ)
OBJECTS_DEBUG = $(patsubst $(SRCDIR)/%.f90,$(OBJDIR)/debug_%.o,$(SOURCES)) $(MODEL_OBJ_DEBUG)

# Phony targets
.PHONY: all debug profile clean distclean test help rnld template models

# Default target
all: $(BINARY)

# Main executable
$(BINARY): $(OBJECTS) | $(BINDIR)
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

# Debug build (model-suffixed)
debug: $(BINDIR)/mass_inf-$(MODEL)-debug

$(BINDIR)/mass_inf-$(MODEL)-debug: $(OBJECTS_DEBUG) | $(BINDIR)
	$(FC) $(FCFLAGS_DEBUG) -o $@ $^ $(LDFLAGS)

# Profile build (model-suffixed)
profile: clean $(BINDIR)/mass_inf-$(MODEL)-prof

$(BINDIR)/mass_inf-$(MODEL)-prof: $(OBJECTS) | $(BINDIR)
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

$(OBJDIR)/debug_grid_config.o: $(OBJDIR)/debug_precision.o
$(MODEL_OBJ_DEBUG): $(OBJDIR)/debug_precision.o $(OBJDIR)/debug_grid_config.o
$(OBJDIR)/debug_pde_stepper.o: $(OBJDIR)/debug_precision.o $(MODEL_OBJ_DEBUG)
$(OBJDIR)/debug_evolve.o: $(OBJDIR)/debug_precision.o $(OBJDIR)/debug_pde_stepper.o $(MODEL_OBJ_DEBUG)
$(OBJDIR)/debug_mass_inf.o: $(OBJDIR)/debug_precision.o $(MODEL_OBJ_DEBUG) $(OBJDIR)/debug_evolve.o $(OBJDIR)/debug_grid_config.o $(OBJDIR)/debug_polint.o $(OBJDIR)/debug_utils.o $(OBJDIR)/debug_AMR.o
$(OBJDIR)/debug_AMR.o: $(OBJDIR)/debug_precision.o $(OBJDIR)/debug_grid_config.o $(OBJDIR)/debug_polint.o $(OBJDIR)/debug_utils.o $(MODEL_OBJ_DEBUG)

# Create directories
$(BINDIR) $(OBJDIR):
	mkdir -p $@

# Clean build artifacts (keep binaries)
clean:
	rm -rf $(OBJDIR) *.mod *.MOD gmon.out

# Full cleanup (removes binaries and all generated files)
distclean: clean
	rm -rf $(BINDIR) timing.log gprof.txt
	rm -rf TESTING

# Run tests
test: all
	@./test/run_tests.sh

# Convenience targets to build specific models with distinct binaries
rnld:
	$(MAKE) --no-print-directory MODEL=rnld all

# Build the following models
models:
	$(MAKE) --no-print-directory clean
	$(MAKE) --no-print-directory rnld
	$(MAKE) --no-print-directory clean

# Help
help:
	@echo "Available targets:"
	@echo "  make                         - Build binary for default model ($(MODEL)) -> bin/mass_inf-$(MODEL)"
	@echo "  make MODEL=<model>           - Build given model (e.g., make MODEL=flat)"
	@echo "  make rnld                    - Build RNLD model -> bin/mass_inf-rnld"
	@echo "  make models                  - Build several models"
	@echo "  make [MODEL=<model>] debug   - Build model with debug flags -> bin/mass_inf-<model>-debug"
	@echo "  make [MODEL=<model>] profile - Build model with gprof instrumentation -> bin/mass_inf-<model>-prof"
	@echo "  make [MODEL=<model>] test    - Run regression tests"
	@echo "  make clean                   - Remove build artifacts"
	@echo "  make distclean               - Remove all generated files"
	@echo "  make help                    - Show this message"
