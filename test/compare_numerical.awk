#!/usr/bin/awk -f
#
# Numerical comparison tool for floating-point data files
# Usage: ./compare_numerical.awk -v rtol=1e-5 -v atol=1e-8 reference.dat test.dat
#
# Compares two files column-by-column using:
#   |test - reference| <= atol + rtol * |reference|
#

BEGIN {
    # Default tolerances if not provided
    if (rtol == "") rtol = 1e-5
    if (atol == "") atol = 1e-8
    
    n_diffs = 0
    total_vals = 0
    file_num = 0
    max_show = 10
    n_shown = 0
}

# Track which file we're reading
FNR == 1 { file_num++ }

# Skip comment lines
/^#/ { next }

# Skip empty lines
NF == 0 { next }

# First file (reference)
file_num == 1 {
    for (i = 1; i <= NF; i++) {
        ref[FNR, i] = $i
    }
    ref_rows = FNR
    ref_cols = NF
    next
}

# Second file (test)
file_num == 2 {
    test_rows = FNR
    test_cols = NF
    
    # Compare this row
    for (i = 1; i <= NF; i++) {
        total_vals++
        test_val = $i
        ref_val = ref[FNR, i]
        
        # Calculate absolute difference
        abs_diff = test_val - ref_val
        if (abs_diff < 0) abs_diff = -abs_diff
        
        # Calculate tolerance threshold
        abs_ref = ref_val
        if (abs_ref < 0) abs_ref = -abs_ref
        threshold = atol + rtol * abs_ref
        
        # Check if within tolerance
        if (abs_diff > threshold) {
            n_diffs++
            col_diffs[i]++
            
            # Track max differences per column
            if (abs_diff > max_abs_diff[i]) {
                max_abs_diff[i] = abs_diff
            }
            
            # Relative difference
            rel_diff = abs_diff / (abs_ref + 1e-100)
            if (rel_diff > max_rel_diff[i]) {
                max_rel_diff[i] = rel_diff
            }
            
            # Store first few for detailed output
            if (n_shown < max_show) {
                diff_row[n_shown] = FNR
                diff_col[n_shown] = i
                diff_ref[n_shown] = ref_val
                diff_test[n_shown] = test_val
                diff_abs[n_shown] = abs_diff
                diff_rel[n_shown] = rel_diff
                n_shown++
            }
        }
    }
}

END {
    # Check for shape mismatch
    if (ref_rows != test_rows || ref_cols != test_cols) {
        printf "Shape mismatch: reference (%d, %d) vs test (%d, %d)\n", \
            ref_rows, ref_cols, test_rows, test_cols
        exit 1
    }
    
    # Report results
    if (n_diffs == 0) {
        printf "✓ All %d rows and %d columns match within tolerance\n", \
            ref_rows, ref_cols
        exit 0
    }
    
    # Differences found
    printf "✗ Numerical differences detected:\n\n"
    printf "  Total differences: %d out of %d values\n\n", n_diffs, total_vals
    
    # Per-column statistics
    printf "  Per-column statistics:\n"
    for (i = 1; i <= ref_cols; i++) {
        if (col_diffs[i] > 0) {
            printf "    col_%d     : %6d diffs, max_abs=%.3e, max_rel=%.3e\n", \
                i, col_diffs[i], max_abs_diff[i], max_rel_diff[i]
        }
    }
    printf "\n"
    
    # Show first few differences
    printf "  First %d differences:\n", n_shown
    printf "    %6s %6s %15s %15s %12s %12s\n", \
        "Row", "Col", "Reference", "Test", "Abs Diff", "Rel Diff"
    
    for (i = 0; i < n_shown; i++) {
        printf "    %6d %6d %15.6e %15.6e %12.3e %12.3e\n", \
            diff_row[i], diff_col[i], diff_ref[i], diff_test[i], \
            diff_abs[i], diff_rel[i]
    }
    
    exit 1
}
