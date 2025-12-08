#!/bin/bash
#
# Regression test runner for mass_inf simulation
#
# This script runs the simulation and compares output against reference data
# using numerical tolerances appropriate for floating-point computation.
#

set -e  # Exit on error

# Activate conda environment if available
if command -v conda &> /dev/null; then
    eval "$(conda shell.bash hook)"
    conda activate numpy 2>/dev/null || true
fi

# Configuration
TEST_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJ_DIR="$(dirname "$TEST_DIR")"
BIN_DIR="$PROJ_DIR/bin"
EXECUTABLE="$BIN_DIR/mass_inf"

REFERENCE_FILE="$TEST_DIR/reference_data.dat"
OUTPUT_FILE="$TEST_DIR/output_data.dat"
DATA_FILE="$PROJ_DIR/data.dat"

# Numerical tolerances
RTOL="1e-5"  # Relative tolerance (0.001%)
ATOL="2e-8"  # Absolute tolerance

# Column names for better reporting (avoid bash special var COLUMNS)
COLUMN_NAMES="u v r phi sigma mass drdv Ricci"

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo "=========================================="
echo "  Mass Inflation - Regression Tests"
echo "=========================================="
echo

# Check if executable exists
if [ ! -f "$EXECUTABLE" ]; then
    echo -e "${RED}Error: Executable not found: $EXECUTABLE${NC}"
    echo "Run 'make release' first to build the code."
    exit 1
fi

# Check if reference data exists
if [ ! -f "$REFERENCE_FILE" ]; then
    echo -e "${YELLOW}Warning: No reference data found.${NC}"
    echo "This appears to be the first run. Creating reference data..."
    echo
    
    # Run simulation
    cd "$PROJ_DIR"
    "$EXECUTABLE" > /dev/null 2>&1
    
    # Copy data.dat to reference
    if [ -f "$DATA_FILE" ]; then
        cp "$DATA_FILE" "$REFERENCE_FILE"
        echo -e "${GREEN}✓ Reference data saved to: $REFERENCE_FILE${NC}"
        echo
        echo "You can now make changes and run this test again to validate."
        exit 0
    else
        echo -e "${RED}Error: Simulation did not produce data.dat${NC}"
        exit 1
    fi
fi

# Run the simulation
echo "Running simulation..."
cd "$PROJ_DIR"
"$EXECUTABLE" > /dev/null 2>&1

# Check if output was created
if [ ! -f "$DATA_FILE" ]; then
    echo -e "${RED}Error: Simulation did not produce data.dat${NC}"
    exit 1
fi

# Copy to test directory for comparison
cp "$DATA_FILE" "$OUTPUT_FILE"

echo -e "${GREEN}✓ Simulation completed${NC}"
echo

# Run numerical comparison
echo "Comparing output with reference data..."
echo "  Relative tolerance: $RTOL"
echo "  Absolute tolerance: $ATOL"
echo

if python3 "$TEST_DIR/compare_numerical.py" \
    "$REFERENCE_FILE" "$OUTPUT_FILE" \
    --rtol "$RTOL" --atol "$ATOL" \
    --columns $COLUMN_NAMES; then
    
    echo
    echo -e "${GREEN}=========================================="
    echo -e "  ALL TESTS PASSED"
    echo -e "==========================================${NC}"
    exit 0
else
    echo
    echo -e "${RED}=========================================="
    echo -e "  TESTS FAILED"
    echo -e "==========================================${NC}"
    echo
    echo "If the changes are intentional, update reference data with:"
    echo "  cp test/output_data.dat test/reference_data.dat"
    exit 1
fi
