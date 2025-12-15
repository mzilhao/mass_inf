FC = gfortran
FCFLAGS = -O3 -march=native -fno-unsafe-math-optimizations
FCFLAGS_DEBUG = -g -Wall -Wextra -Wno-unused-dummy-argument -fcheck=all \
                -ffpe-trap=invalid,zero,overflow -Wno-maybe-uninitialized
FCFLAGS_PROFILE = -O2 -pg
LDFLAGS =

# Directories
BINDIR = bin
OBJDIR = obj

# Source files and objects
SOURCES = $(wildcard *.f90)
OBJECTS = $(patsubst %.f90,$(OBJDIR)/%.o,$(SOURCES))
OBJECTS_DEBUG = $(patsubst %.f90,$(OBJDIR)/debug_%.o,$(SOURCES))

# Phony targets
.PHONY: all debug profile clean distclean test test-save-reference help

# Default target
all: $(BINDIR)/mass_inf

# Main executable
$(BINDIR)/mass_inf: $(OBJECTS) | $(BINDIR)
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

# Compile release objects
$(OBJDIR)/%.o: %.f90 | $(OBJDIR)
	$(FC) $(FCFLAGS) -J$(OBJDIR) -c -o $@ $<

# Compile debug objects
$(OBJDIR)/debug_%.o: %.f90 | $(OBJDIR)
	$(FC) $(FCFLAGS_DEBUG) -J$(OBJDIR) -c -o $@ $<

# Module dependencies
$(OBJDIR)/functions.o: $(OBJDIR)/simulation_config.o
$(OBJDIR)/pde_stepper.o: $(OBJDIR)/functions.o
$(OBJDIR)/evolve.o: $(OBJDIR)/pde_stepper.o $(OBJDIR)/functions.o
$(OBJDIR)/mass_inf.o: $(OBJDIR)/functions.o $(OBJDIR)/evolve.o $(OBJDIR)/simulation_config.o $(OBJDIR)/polint.o $(OBJDIR)/utils.o

$(OBJDIR)/debug_functions.o: $(OBJDIR)/debug_simulation_config.o
$(OBJDIR)/debug_pde_stepper.o: $(OBJDIR)/debug_functions.o
$(OBJDIR)/debug_evolve.o: $(OBJDIR)/debug_pde_stepper.o $(OBJDIR)/debug_functions.o
$(OBJDIR)/debug_mass_inf.o: $(OBJDIR)/debug_functions.o $(OBJDIR)/debug_evolve.o $(OBJDIR)/debug_simulation_config.o $(OBJDIR)/debug_polint.o $(OBJDIR)/debug_utils.o

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
	@echo "  make clean     - Remove build artifacts"
	@echo "  make distclean - Remove all generated files"
	@echo "  make help      - Show this message"
