"""
Data Validation Script for Meta-Analysis
Run this before your main analysis to check data format
"""

import pandas as pd
import numpy as np

def validate_data(csv_path):
    """
    Validate the input CSV file for meta-analysis
    
    Parameters:
    -----------
    csv_path : str
        Path to the CSV file
    """
    print("="*70)
    print("DATA VALIDATION FOR META-ANALYSIS")
    print("="*70)
    
    # Required columns
    required_cols = [
        'Study_ID', 'Setting', 'Sensitivity', 'Specificity', 
        'AUC', 'n_diseased', 'n_non_diseased', 'Device', 'Reference_Test'
    ]
    
    try:
        # Load data
        df = pd.read_csv(csv_path)
        print(f"\n✓ Successfully loaded file: {csv_path}")
        print(f"  Rows: {len(df)}")
        print(f"  Columns: {len(df.columns)}")
        
    except FileNotFoundError:
        print(f"\n✗ ERROR: File not found: {csv_path}")
        print("  Please check the file path and name.")
        return False
    except Exception as e:
        print(f"\n✗ ERROR: Could not read file: {str(e)}")
        return False
    
    # Check columns
    print("\n" + "-"*70)
    print("CHECKING COLUMNS")
    print("-"*70)
    
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        print(f"\n✗ Missing required columns: {missing_cols}")
        print(f"\n  Current columns: {list(df.columns)}")
        print(f"  Required columns: {required_cols}")
        return False
    else:
        print("✓ All required columns present")
    
    # Check data types and values
    print("\n" + "-"*70)
    print("CHECKING DATA QUALITY")
    print("-"*70)
    
    issues = []
    
    # Check for missing values
    missing_counts = df[required_cols].isnull().sum()
    if missing_counts.any():
        print("\n⚠ WARNING: Missing values detected:")
        for col, count in missing_counts[missing_counts > 0].items():
            print(f"  {col}: {count} missing values")
            issues.append(f"Missing values in {col}")
    else:
        print("✓ No missing values")
    
    # Check Sensitivity range (should be 0-1)
    if 'Sensitivity' in df.columns:
        sens_min, sens_max = df['Sensitivity'].min(), df['Sensitivity'].max()
        if sens_min < 0 or sens_max > 1:
            if sens_max > 1 and sens_max <= 100:
                print("\n⚠ WARNING: Sensitivity appears to be in percentage format (0-100)")
                print("  Please convert to proportions (0-1)")
                issues.append("Sensitivity in wrong format")
            else:
                print(f"\n✗ ERROR: Sensitivity out of range: {sens_min:.3f} - {sens_max:.3f}")
                issues.append("Sensitivity out of range")
        else:
            print(f"✓ Sensitivity range valid: {sens_min:.3f} - {sens_max:.3f}")
    
    # Check Specificity range (should be 0-1)
    if 'Specificity' in df.columns:
        spec_min, spec_max = df['Specificity'].min(), df['Specificity'].max()
        if spec_min < 0 or spec_max > 1:
            if spec_max > 1 and spec_max <= 100:
                print("\n⚠ WARNING: Specificity appears to be in percentage format (0-100)")
                print("  Please convert to proportions (0-1)")
                issues.append("Specificity in wrong format")
            else:
                print(f"\n✗ ERROR: Specificity out of range: {spec_min:.3f} - {spec_max:.3f}")
                issues.append("Specificity out of range")
        else:
            print(f"✓ Specificity range valid: {spec_min:.3f} - {spec_max:.3f}")
    
    # Check AUC range (should be 0.5-1)
    if 'AUC' in df.columns:
        auc_min, auc_max = df['AUC'].min(), df['AUC'].max()
        if auc_min < 0.5 or auc_max > 1:
            print(f"\n⚠ WARNING: AUC out of typical range: {auc_min:.3f} - {auc_max:.3f}")
            print("  AUC should typically be between 0.5 and 1.0")
            issues.append("AUC unusual range")
        else:
            print(f"✓ AUC range valid: {auc_min:.3f} - {auc_max:.3f}")
    
    # Check sample sizes
    if 'n_diseased' in df.columns:
        if (df['n_diseased'] < 1).any():
            print("\n✗ ERROR: n_diseased contains values < 1")
            issues.append("Invalid n_diseased")
        else:
            print(f"✓ n_diseased valid (range: {df['n_diseased'].min()}-{df['n_diseased'].max()})")
    
    if 'n_non_diseased' in df.columns:
        if (df['n_non_diseased'] < 1).any():
            print("\n✗ ERROR: n_non_diseased contains values < 1")
            issues.append("Invalid n_non_diseased")
        else:
            print(f"✓ n_non_diseased valid (range: {df['n_non_diseased'].min()}-{df['n_non_diseased'].max()})")
    
    # Check Setting values
    if 'Setting' in df.columns:
        settings = df['Setting'].unique()
        print(f"\n✓ Settings found: {list(settings)}")
        print(f"  Counts: {df['Setting'].value_counts().to_dict()}")
        
        if len(settings) < 2:
            print("  ⚠ Note: Only one setting - subgroup analysis won't be meaningful")
    
    # Check for duplicate Study_IDs
    if 'Study_ID' in df.columns:
        duplicates = df['Study_ID'].duplicated()
        if duplicates.any():
            dup_ids = df[duplicates]['Study_ID'].tolist()
            print(f"\n✗ ERROR: Duplicate Study_IDs found: {dup_ids}")
            issues.append("Duplicate Study_IDs")
        else:
            print("✓ All Study_IDs unique")
    
    # Summary
    print("\n" + "="*70)
    print("VALIDATION SUMMARY")
    print("="*70)
    
    if not issues:
        print("\n✓✓✓ DATA IS VALID AND READY FOR ANALYSIS! ✓✓✓")
        print("\nYou can now run the meta-analysis with:")
        print("  python run_analysis.py")
        print("\nOr in Python:")
        print("  from diagnostic_meta_analysis import DiagnosticMetaAnalysis")
        print("  meta = DiagnosticMetaAnalysis('niri_meta.csv')")
        print("  results = meta.run_complete_analysis()")
        return True
    else:
        print(f"\n⚠ ISSUES FOUND: {len(issues)}")
        for i, issue in enumerate(issues, 1):
            print(f"  {i}. {issue}")
        print("\nPlease fix these issues before running the analysis.")
        return False

# Additional helper function to preview data
def preview_data(csv_path, n_rows=5):
    """Preview the first n rows of the dataset"""
    df = pd.read_csv(csv_path)
    print("\n" + "="*70)
    print(f"DATA PREVIEW (first {n_rows} rows)")
    print("="*70)
    print(df.head(n_rows).to_string())
    print("\n" + "="*70)
    print("BASIC STATISTICS")
    print("="*70)
    numeric_cols = ['Sensitivity', 'Specificity', 'AUC', 'n_diseased', 'n_non_diseased']
    print(df[numeric_cols].describe())


if __name__ == "__main__":
    import sys
    
    # Get filename from command line or use default
    if len(sys.argv) > 1:
        csv_file = sys.argv[1]
    else:
        csv_file = 'niri_meta.csv'
    
    # Run validation
    is_valid = validate_data(csv_file)
    
    # Show preview if valid
    if is_valid:
        try:
            preview_data(csv_file)
        except:
            pass
    
    print("\n" + "="*70)
