"""
Simple Example Usage for Diagnostic Meta-Analysis
Run this script after installing required packages
"""

from diagnostic_meta_analysis import DiagnosticMetaAnalysis

# Initialize with your CSV file
meta = DiagnosticMetaAnalysis('niri_meta.csv')

# Run the complete analysis
# This will:
# 1. Prepare the data (calculate TP, FN, FP, TN)
# 2. Perform overall bivariate meta-analysis
# 3. Perform subgroup analysis (in vivo vs in vitro)
# 4. Create forest plots with subgroups
# 5. Create SROC curve
# 6. Create coupled forest plot
# 7. Generate summary table

results = meta.run_complete_analysis()

# Access results
print("\n" + "="*60)
print("YOUR RESULTS:")
print("="*60)
print(f"Pooled Sensitivity: {results['overall_sensitivity']:.3f} ({results['overall_sensitivity']*100:.1f}%)")
print(f"Pooled Specificity: {results['overall_specificity']:.3f} ({results['overall_specificity']*100:.1f}%)")

# View subgroup results
print("\nSubgroup Analysis:")
for setting, data in results['subgroup_results'].items():
    if data:
        print(f"\n{setting}:")
        print(f"  Sensitivity: {data['sensitivity']:.3f} (95% CI: {data['sens_ci'][0]:.3f}-{data['sens_ci'][1]:.3f})")
        print(f"  Specificity: {data['specificity']:.3f} (95% CI: {data['spec_ci'][0]:.3f}-{data['spec_ci'][1]:.3f})")

print("\n✓ All figures saved to /mnt/user-data/outputs/")
