#!/bin/bash
#
# Regression test runner for mass_inf simulation
#
# This script runs the simulation and compares output against reference data
# using numerical tolerances appropriate for floating-point computation.

set -u  # Treat unset vars as error

# Activate conda environment if available
if command -v conda &> /dev/null; then
    eval "$(conda shell.bash hook)"
    conda activate numpy 2>/dev/null || true
fi

TEST_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJ_DIR="$(dirname "$TEST_DIR")"
BIN_DIR="$PROJ_DIR/bin"

# Model selection (must be set by Makefile or environment)
if [ -z "${MODEL:-}" ]; then
    echo -e "${RED}Error: MODEL variable not set. Please run with 'make test MODEL=<model>' or set MODEL in the environment.${NC}"
    exit 1
fi
MODEL="$MODEL"

# Prefer model-suffixed binary; fall back to unsuffixed if present
EXECUTABLE="$BIN_DIR/mass_inf-$MODEL"
if [ ! -x "$EXECUTABLE" ] && [ -x "$BIN_DIR/mass_inf" ]; then
    EXECUTABLE="$BIN_DIR/mass_inf"
fi

# Output root (per model to avoid clashes)
OUT_ROOT="$PROJ_DIR/TESTING/$MODEL"

# Numerical tolerances
RTOL="1e-5"  # Relative tolerance (0.001%)
ATOL="2e-8"  # Absolute tolerance

# Column names for reporting (match current outputs)
COLUMNS_FIELDS="u r phi sigma"
COLUMNS_DERIVS="u dr_du dr_dv"
COLUMNS_DIAG="u mass ricci"

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Discover test cases: require model-specific folder test/<model>
CASE_ROOT="$TEST_DIR/$MODEL"
if [ ! -d "$CASE_ROOT" ]; then
    echo -e "${YELLOW}No test directory for model '$MODEL': $CASE_ROOT${NC}"
    echo -e "${YELLOW}Add configs and references under test/$MODEL/ or choose a model with tests.${NC}"
    exit 0
fi

# Print banner after resolving paths
echo "=========================================="
echo "  Mass Inflation - Regression Tests"
echo "=========================================="
echo "Model:       $MODEL"
echo "Executable:  $EXECUTABLE"
echo "Cases root:  $CASE_ROOT"
echo "Output root: $OUT_ROOT"
echo

# Check if executable exists
if [ ! -f "$EXECUTABLE" ]; then
    echo -e "${RED}Error: Executable not found: $EXECUTABLE${NC}"
    echo "Run 'make release' first to build the code."
    exit 1
fi

# Helper to compare two files with column labels
compare_file() {
    local ref_file="$1"
    local out_file="$2"
    local columns="$3"
    python3 "$TEST_DIR/compare_numerical.py" "$ref_file" "$out_file" --rtol "$RTOL" --atol "$ATOL" --columns $columns
}

shopt -s nullglob
cases=()
for cfg in "$CASE_ROOT"/*.nml; do
    cases+=("$cfg")
done
shopt -u nullglob

if [ ${#cases[@]} -eq 0 ]; then
    echo -e "${YELLOW}No test configurations (*.nml) found in $CASE_ROOT${NC}"
    exit 0
fi

overall_status=0

for cfg in "${cases[@]}"; do
    case_name="$(basename "$cfg" .nml)"

    # Reference directory sits alongside the config (either model-specific or root)
    ref_dir="$(dirname "$cfg")/$case_name"
    out_dir="$OUT_ROOT/$case_name"

    echo "------------------------------------------"
    echo "Test case: $case_name"
    echo "Config:    $cfg"
    echo "Ref dir:   $ref_dir"
    echo "Out dir:   $out_dir"

    # Prepare output root
    mkdir -p "$OUT_ROOT"
    # Clean previous output
    rm -rf "$out_dir"

    echo "Running simulation..."
    # Run from the per-model OUT_ROOT so outputs land under that tree
    if ! (cd "$OUT_ROOT" && "$EXECUTABLE" "$cfg" > /dev/null 2>&1); then
        echo -e "${RED}Execution failed for $case_name${NC}"
        overall_status=1
        continue
    fi

    # Verify outputs
    missing=0
    for f in fields.dat derivatives.dat diagnostics.dat; do
        if [ ! -f "$out_dir/$f" ]; then
            echo -e "${RED}Missing output file: $out_dir/$f${NC}"
            missing=1
        fi
    done
    if [ "$missing" -ne 0 ]; then
        overall_status=1
        continue
    fi

    # Ensure reference directory exists with expected files
    if [ ! -d "$ref_dir" ]; then
        echo -e "${RED}Missing reference directory: $ref_dir${NC}"
        overall_status=1
        continue
    fi
    missing_ref=0
    for f in fields.dat derivatives.dat diagnostics.dat; do
        if [ ! -f "$ref_dir/$f" ]; then
            echo -e "${RED}Missing reference file: $ref_dir/$f${NC}"
            missing_ref=1
        fi
    done
    if [ "$missing_ref" -ne 0 ]; then
        overall_status=1
        continue
    fi

    echo "Comparing outputs..."
    status_case=0
    if ! compare_file "$ref_dir/fields.dat" "$out_dir/fields.dat" "$COLUMNS_FIELDS"; then
        status_case=1
    fi
    if ! compare_file "$ref_dir/derivatives.dat" "$out_dir/derivatives.dat" "$COLUMNS_DERIVS"; then
        status_case=1
    fi
    if ! compare_file "$ref_dir/diagnostics.dat" "$out_dir/diagnostics.dat" "$COLUMNS_DIAG"; then
        status_case=1
    fi

    if [ "$status_case" -eq 0 ]; then
        echo -e "${GREEN}✓ Test passed: $case_name${NC}"
    else
        echo -e "${RED}✗ Test failed: $case_name${NC}"
        overall_status=1
    fi
done

echo
if [ "$overall_status" -eq 0 ]; then
    echo -e "${GREEN}=========================================="
    echo -e "  ALL TESTS PASSED"
    echo -e "==========================================${NC}"
    exit 0
else
    echo -e "${RED}=========================================="
    echo -e "  TESTS FAILED"
    echo -e "==========================================${NC}"
    echo "Update reference outputs if changes are expected."
    exit 1
fi
