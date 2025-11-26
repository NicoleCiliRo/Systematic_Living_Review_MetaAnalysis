"""
CSV File Diagnostic Tool
Identifies formatting issues in your meta-analysis CSV file
"""

import csv

def diagnose_csv(filename):
    """
    Diagnose issues in CSV file
    """
    print("="*70)
    print("CSV FILE DIAGNOSTIC TOOL")
    print("="*70)
    print(f"\nAnalyzing: {filename}\n")
    
    try:
        with open(filename, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        
        print(f"✓ File opened successfully")
        print(f"✓ Total lines: {len(lines)}\n")
        
        # Check each line
        print("-"*70)
        print("LINE-BY-LINE ANALYSIS:")
        print("-"*70)
        
        expected_cols = None
        issues_found = False
        
        for i, line in enumerate(lines, 1):
            # Count fields (simple comma count + 1)
            # This is rough but helps identify issues
            line_stripped = line.strip()
            
            # Skip empty lines
            if not line_stripped:
                print(f"Line {i}: EMPTY LINE (should be removed)")
                issues_found = True
                continue
            
            # Count commas (fields = commas + 1, approximately)
            comma_count = line_stripped.count(',')
            field_count = comma_count + 1
            
            # First line sets expected column count
            if i == 1:
                expected_cols = field_count
                print(f"Line {i} (HEADER): {field_count} fields")
                print(f"  Preview: {line_stripped[:80]}...")
                print(f"  ✓ This sets the expected column count: {expected_cols}")
            else:
                status = "✓" if field_count == expected_cols else "✗ ERROR"
                print(f"Line {i}: {field_count} fields {status}")
                
                if field_count != expected_cols:
                    print(f"  Expected: {expected_cols} fields")
                    print(f"  Found: {field_count} fields")
                    print(f"  Preview: {line_stripped[:100]}...")
                    issues_found = True
                    
                    # Try to identify the issue
                    if field_count > expected_cols:
                        print(f"  ⚠ EXTRA FIELD(S) - Check for:")
                        print(f"     • Unquoted commas in text fields")
                        print(f"     • Extra columns accidentally added")
                    else:
                        print(f"  ⚠ MISSING FIELD(S) - Check for:")
                        print(f"     • Missing values (use empty quotes \"\")")
                        print(f"     • Deleted columns")
                    print()
        
        # Summary
        print("="*70)
        if not issues_found:
            print("✓✓✓ NO ISSUES FOUND! Your CSV format looks good.")
            print("\nTry running validation again:")
            print("  python validate_data.py niri_meta.csv")
        else:
            print("⚠ ISSUES FOUND - See details above")
            print("\nCommon fixes:")
            print("  1. Remove empty lines")
            print("  2. Put quotes around text with commas: \"Device A, Model X\"")
            print("  3. Ensure all rows have same number of columns")
            print("  4. Check for hidden characters or extra spaces")
        print("="*70)
        
    except FileNotFoundError:
        print(f"✗ ERROR: File '{filename}' not found")
        print("  Make sure the file is in the current directory")
    except Exception as e:
        print(f"✗ ERROR: {str(e)}")

def show_csv_structure(filename):
    """
    Show the actual parsed structure using Python's csv module
    """
    print("\n" + "="*70)
    print("DETAILED CSV PARSING:")
    print("="*70)
    
    try:
        with open(filename, 'r', encoding='utf-8') as f:
            reader = csv.reader(f)
            for i, row in enumerate(reader, 1):
                print(f"\nLine {i}: {len(row)} fields")
                for j, field in enumerate(row, 1):
                    print(f"  Field {j}: '{field}'")
                if i >= 5:  # Only show first 5 rows
                    print("\n  ... (showing first 5 rows only)")
                    break
    except Exception as e:
        print(f"Could not parse CSV: {str(e)}")

def fix_common_issues(filename, output_filename=None):
    """
    Attempt to auto-fix common issues
    """
    if output_filename is None:
        output_filename = filename.replace('.csv', '_fixed.csv')
    
    print("\n" + "="*70)
    print("ATTEMPTING AUTO-FIX:")
    print("="*70)
    
    try:
        with open(filename, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        
        # Remove empty lines
        lines = [line for line in lines if line.strip()]
        
        # Remove any BOM or hidden characters
        lines = [line.replace('\ufeff', '') for line in lines]
        
        with open(output_filename, 'w', encoding='utf-8', newline='') as f:
            f.writelines(lines)
        
        print(f"✓ Created cleaned file: {output_filename}")
        print(f"✓ Removed empty lines")
        print(f"✓ Removed hidden characters")
        print(f"\nNow try validating the fixed file:")
        print(f"  python validate_data.py {output_filename}")
        
    except Exception as e:
        print(f"Could not auto-fix: {str(e)}")

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        filename = 'niri_meta.csv'
    
    print("\n")
    diagnose_csv(filename)
    show_csv_structure(filename)
    
    # Offer to create a fixed version
    print("\n" + "="*70)
    response = input("\nWould you like to create an auto-fixed version? (y/n): ")
    if response.lower() in ['y', 'yes']:
        fix_common_issues(filename)
    
    print("\n")
