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
.PHONY: all debug profile clean distclean test help rnld flat template models rnld-profile flat-profile rnld-debug flat-debug

# Default target
all: $(BINARY)

# Main executable
$(BINARY): $(OBJECTS) | $(BINDIR)
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

# Debug build (model-suffixed)
debug: $(BINDIR)/mass_inf-$(MODEL)-debug

$(BINDIR)/mass_inf-$(MODEL)-debug: $(OBJECTS_DEBUG) | $(BINDIR)
	$(FC) $(FCFLAGS_DEBUG) -o $@ $^ $(LDFLAGS)

# Debug convenience targets per model (name: <model>-debug)
rnld-debug:
	$(MAKE) --no-print-directory MODEL=rnld debug

flat-debug:
	$(MAKE) --no-print-directory MODEL=flat debug

# Profile build (model-suffixed)
profile: clean $(BINDIR)/mass_inf-$(MODEL)-prof

$(BINDIR)/mass_inf-$(MODEL)-prof: $(OBJECTS) | $(BINDIR)
	$(FC) $(FCFLAGS_PROFILE) -o $@ $^ -pg

# Profile convenience targets per model (name: <model>-profile)
rnld-profile:
	$(MAKE) --no-print-directory MODEL=rnld profile

flat-profile:
	$(MAKE) --no-print-directory MODEL=flat profile

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

# Clean build artifacts
clean:
	rm -rf $(OBJDIR) $(BINDIR) *.mod *.MOD gmon.out

# Full cleanup
distclean: clean
	rm -f timing.log gprof.txt
	rm -rf TESTING

# Run tests
test: all
	@./test/run_tests.sh

# Convenience targets to build specific models with distinct binaries
rnld:
	$(MAKE) --no-print-directory MODEL=rnld BIN_SUFFIX=-rnld all

flat:
	$(MAKE) --no-print-directory MODEL=flat BIN_SUFFIX=-flat all

template:
	$(MAKE) --no-print-directory MODEL=template BIN_SUFFIX=-template all

# Build supported models (template is excluded; it's a scaffolding file)
models:
	$(MAKE) --no-print-directory clean
	$(MAKE) --no-print-directory rnld
	$(MAKE) --no-print-directory clean
	$(MAKE) --no-print-directory flat

# Help
help:
	@echo "Available targets:"
	@echo "  make               - Build optimized binary for default model ($(MODEL)) -> bin/mass_inf-$(MODEL)"
	@echo "  make rnld          - Build RNLD model -> bin/mass_inf-rnld"
	@echo "  make flat          - Build flat model -> bin/mass_inf-flat"
	@echo "  make template      - Build template model -> bin/mass_inf-template (scaffolding only)"
	@echo "  make models        - Build supported models (rnld, flat)"
	@echo "  make debug         - Build with debug flags and checks -> bin/mass_inf-$(MODEL)-debug"
	@echo "  make profile       - Build with gprof instrumentation -> bin/mass_inf-$(MODEL)-prof"
	@echo "  make rnld-debug    - Debug RNLD model -> bin/mass_inf-rnld-debug"
	@echo "  make flat-debug    - Debug flat model -> bin/mass_inf-flat-debug"
	@echo "  make rnld-profile  - Profile RNLD model -> bin/mass_inf-rnld-prof"
	@echo "  make flat-profile  - Profile flat model -> bin/mass_inf-flat-prof"
	@echo "  make test          - Run regression tests"
	@echo "  make clean         - Remove build artifacts"
	@echo "  make distclean     - Remove all generated files"
	@echo "  make help          - Show this message"
