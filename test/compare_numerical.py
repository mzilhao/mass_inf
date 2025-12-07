#!/usr/bin/env python3
"""
Numerical comparison tool for regression testing of Fortran simulation output.

Compares two data files containing multi-column floating-point output,
using configurable absolute and relative tolerances to account for
numerical precision limitations.
"""

import sys
import argparse
import numpy as np


def parse_data_file(filename):
    """
    Parse a data file, skipping comment lines starting with '#'.
    
    Returns:
        numpy.ndarray: Data array with shape (n_rows, n_columns)
    """
    try:
        data = np.loadtxt(filename, comments='#')
        return data
    except Exception as e:
        print(f"Error reading file {filename}: {e}")
        sys.exit(1)


def compare_arrays(reference, test, rtol=1e-5, atol=1e-8, column_names=None):
    """
    Compare two numpy arrays with relative and absolute tolerance.
    
    Args:
        reference: Reference (golden) data array
        test: Test data array to validate
        rtol: Relative tolerance (default: 1e-5, i.e., 0.001%)
        atol: Absolute tolerance (default: 1e-8)
        column_names: Optional list of column names for better reporting
        
    Returns:
        tuple: (passed, report_string)
    """
    # Check shape consistency
    if reference.shape != test.shape:
        return False, f"Shape mismatch: reference {reference.shape} vs test {test.shape}"
    
    n_rows, n_cols = reference.shape
    
    # Use numpy's allclose for element-wise comparison
    # |test - reference| <= atol + rtol * |reference|
    all_close = np.allclose(test, reference, rtol=rtol, atol=atol)
    
    if all_close:
        return True, f"✓ All {n_rows} rows and {n_cols} columns match within tolerance"
    
    # If not all close, provide detailed diagnostics
    report_lines = ["✗ Numerical differences detected:\n"]
    
    # Find which elements differ
    diff_mask = ~np.isclose(test, reference, rtol=rtol, atol=atol)
    diff_positions = np.argwhere(diff_mask)
    
    n_diffs = len(diff_positions)
    report_lines.append(f"  Total differences: {n_diffs} out of {n_rows * n_cols} values\n")
    
    # Show summary statistics for each column
    report_lines.append("  Per-column statistics:")
    for col in range(n_cols):
        col_diffs = diff_mask[:, col].sum()
        if col_diffs > 0:
            col_name = column_names[col] if column_names else f"col_{col}"
            max_abs_diff = np.max(np.abs(test[:, col] - reference[:, col]))
            max_rel_diff = np.max(np.abs((test[:, col] - reference[:, col]) / 
                                        (reference[:, col] + 1e-100)))  # avoid div by zero
            report_lines.append(f"    {col_name:10s}: {col_diffs:6d} diffs, "
                              f"max_abs={max_abs_diff:.3e}, max_rel={max_rel_diff:.3e}")
    
    report_lines.append("")
    
    # Show first few differences in detail
    n_show = min(10, n_diffs)
    report_lines.append(f"  First {n_show} differences:")
    report_lines.append(f"    {'Row':>6s} {'Col':>6s} {'Reference':>15s} {'Test':>15s} {'Abs Diff':>12s} {'Rel Diff':>12s}")
    
    for i in range(n_show):
        row, col = diff_positions[i]
        ref_val = reference[row, col]
        test_val = test[row, col]
        abs_diff = abs(test_val - ref_val)
        rel_diff = abs_diff / (abs(ref_val) + 1e-100)
        
        col_name = column_names[col] if column_names else f"{col}"
        report_lines.append(f"    {row:6d} {col_name:>6s} {ref_val:15.6e} {test_val:15.6e} "
                          f"{abs_diff:12.3e} {rel_diff:12.3e}")
    
    return False, "\n".join(report_lines)


def main():
    parser = argparse.ArgumentParser(
        description="Compare numerical data files for regression testing",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s reference.dat test.dat
  %(prog)s reference.dat test.dat --rtol 1e-6 --atol 1e-10
  %(prog)s reference.dat test.dat --columns u v r phi sigma mass drdv Ricci
        """
    )
    
    parser.add_argument('reference', help='Reference (golden) data file')
    parser.add_argument('test', help='Test data file to validate')
    parser.add_argument('--rtol', type=float, default=1e-5,
                       help='Relative tolerance (default: 1e-5)')
    parser.add_argument('--atol', type=float, default=1e-8,
                       help='Absolute tolerance (default: 1e-8)')
    parser.add_argument('--columns', nargs='+', default=None,
                       help='Optional column names for better reporting')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Verbose output')
    
    args = parser.parse_args()
    
    if args.verbose:
        print(f"Comparing files:")
        print(f"  Reference: {args.reference}")
        print(f"  Test:      {args.test}")
        print(f"  Tolerances: rtol={args.rtol}, atol={args.atol}")
        print()
    
    # Load data files
    reference_data = parse_data_file(args.reference)
    test_data = parse_data_file(args.test)
    
    # Compare
    passed, report = compare_arrays(reference_data, test_data, 
                                    rtol=args.rtol, atol=args.atol,
                                    column_names=args.columns)
    
    print(report)
    
    # Exit with appropriate code
    sys.exit(0 if passed else 1)


if __name__ == '__main__':
    main()
