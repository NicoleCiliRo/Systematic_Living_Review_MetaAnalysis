"""
Diagnostic Test Accuracy Meta-Analysis for Living Systematic Review
Using Python with rpy2 to leverage R's mada and meta packages
Author: Created for NIRI diagnostic accuracy meta-analysis
"""

import pandas as pd
import numpy as np
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
import warnings
warnings.filterwarnings('ignore')

# Activate pandas conversion
pandas2ri.activate()

class DiagnosticMetaAnalysis:
    """
    Class for performing diagnostic test accuracy meta-analysis
    with forest plots, SROC curves, and pooled estimates
    """
    
    def __init__(self, csv_path):
        """
        Initialize with CSV file path
        
        Parameters:
        -----------
        csv_path : str
            Path to CSV file with columns: Study_ID, Setting, Sensitivity, 
            Specificity, AUC, n_diseased, n_non_diseased, Device, Reference_Test
        """
        print("="*70)
        print("DIAGNOSTIC TEST ACCURACY META-ANALYSIS")
        print("="*70)
        
        self.df = pd.read_csv(csv_path)
        print(f"\nLoaded {len(self.df)} studies from {csv_path}")
        self.setup_r_packages()
        
    def setup_r_packages(self):
        """Install and load required R packages"""
        print("\n" + "-"*70)
        print("Setting up R packages...")
        print("-"*70)
        
        # R code to install packages if needed
        ro.r('''
            install_if_missing <- function(pkg) {
                if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
                    cat(paste("Installing", pkg, "...\n"))
                    install.packages(pkg, repos='http://cran.rstudio.com/', quiet = TRUE)
                    library(pkg, character.only = TRUE, quietly = TRUE)
                }
            }
            
            install_if_missing("meta")
            install_if_missing("mada")
            install_if_missing("metafor")
            install_if_missing("ggplot2")
            install_if_missing("dplyr")
        ''')
        
        # Import packages
        self.meta = importr('meta')
        self.mada = importr('mada')
        self.metafor = importr('metafor')
        self.grdevices = importr('grDevices')
        
        print("✓ R packages loaded successfully!\n")
        
    def prepare_data(self):
        """
        Prepare data for meta-analysis
        Calculate TP, FN, FP, TN from sensitivity, specificity, and sample sizes
        """
        print("-"*70)
        print("PREPARING DATA")
        print("-"*70)
        
        # Display original data
        print("\nOriginal data columns:", list(self.df.columns))
        print(f"Number of studies: {len(self.df)}")
        
        # Calculate counts from proportions
        self.df['TP'] = (self.df['Sensitivity'] * self.df['n_diseased']).round().astype(int)
        self.df['FN'] = (self.df['n_diseased'] - self.df['TP']).astype(int)
        self.df['TN'] = (self.df['Specificity'] * self.df['n_non_diseased']).round().astype(int)
        self.df['FP'] = (self.df['n_non_diseased'] - self.df['TN']).astype(int)
        
        # Display breakdown by setting
        print("\nStudies by setting:")
        setting_counts = self.df['Setting'].value_counts()
        for setting, count in setting_counts.items():
            print(f"  {setting}: {count} studies")
        
        # Display summary statistics
        print("\nSummary statistics:")
        print(f"  Sensitivity range: {self.df['Sensitivity'].min():.3f} - {self.df['Sensitivity'].max():.3f}")
        print(f"  Specificity range: {self.df['Specificity'].min():.3f} - {self.df['Specificity'].max():.3f}")
        print(f"  Mean AUC: {self.df['AUC'].mean():.3f}")
        
        return self.df
    
    def perform_overall_meta_analysis(self):
        """
        Perform overall bivariate meta-analysis using REITSMA model
        """
        print("\n" + "="*70)
        print("PERFORMING OVERALL BIVARIATE META-ANALYSIS (REITSMA MODEL)")
        print("="*70)
        
        # Transfer data to R
        r_df = pandas2ri.py2rpy(self.df)
        ro.globalenv['meta_data'] = r_df
        
        # Perform bivariate meta-analysis
        ro.r('''
            # Bivariate meta-analysis using Reitsma model
            fit_reitsma <- reitsma(meta_data, 
                                   correction = 0.5,
                                   correction.control = "all")
            
            # Summary of results
            summary_results <- summary(fit_reitsma)
            
            # Extract pooled estimates (on logit scale)
            pooled_sens_logit <- summary_results$coefficients[1, 1]
            pooled_spec_logit <- summary_results$coefficients[2, 1]
            
            # Standard errors
            se_sens <- summary_results$coefficients[1, 2]
            se_spec <- summary_results$coefficients[2, 2]
            
            # Transform from logit scale to probability
            pooled_sens_trans <- plogis(pooled_sens_logit)
            pooled_spec_trans <- plogis(pooled_spec_logit)
            
            # Calculate 95% CI
            sens_ci_lower <- plogis(pooled_sens_logit - 1.96 * se_sens)
            sens_ci_upper <- plogis(pooled_sens_logit + 1.96 * se_sens)
            spec_ci_lower <- plogis(pooled_spec_logit - 1.96 * se_spec)
            spec_ci_upper <- plogis(pooled_spec_logit + 1.96 * se_spec)
        ''')
        
        # Get results
        pooled_sens = ro.r('pooled_sens_trans')[0]
        pooled_spec = ro.r('pooled_spec_trans')[0]
        sens_ci = (ro.r('sens_ci_lower')[0], ro.r('sens_ci_upper')[0])
        spec_ci = (ro.r('spec_ci_lower')[0], ro.r('spec_ci_upper')[0])
        
        print(f"\n{'POOLED ESTIMATES (Overall)':^70}")
        print("-"*70)
        print(f"Pooled Sensitivity: {pooled_sens:.3f} (95% CI: {sens_ci[0]:.3f} - {sens_ci[1]:.3f})")
        print(f"                    {pooled_sens*100:.1f}% (95% CI: {sens_ci[0]*100:.1f}% - {sens_ci[1]*100:.1f}%)")
        print()
        print(f"Pooled Specificity: {pooled_spec:.3f} (95% CI: {spec_ci[0]:.3f} - {spec_ci[1]:.3f})")
        print(f"                    {pooled_spec*100:.1f}% (95% CI: {spec_ci[0]*100:.1f}% - {spec_ci[1]*100:.1f}%)")
        print("-"*70)
        
        return pooled_sens, pooled_spec
    
    def perform_subgroup_meta_analysis(self):
        """
        Perform subgroup meta-analysis for in vivo and in vitro studies
        """
        print("\n" + "="*70)
        print("PERFORMING SUBGROUP META-ANALYSIS")
        print("="*70)
        
        settings = self.df['Setting'].unique()
        results = {}
        
        for setting in settings:
            print(f"\n{setting.upper()} STUDIES:")
            print("-"*70)
            
            subset = self.df[self.df['Setting'] == setting]
            print(f"Number of studies: {len(subset)}")
            
            # Transfer subset to R
            r_subset = pandas2ri.py2rpy(subset)
            ro.globalenv['subset_data'] = r_subset
            
            try:
                # Perform bivariate meta-analysis on subset
                ro.r('''
                    fit_subset <- reitsma(subset_data, 
                                         correction = 0.5,
                                         correction.control = "all")
                    summary_subset <- summary(fit_subset)
                    
                    pooled_sens_logit <- summary_subset$coefficients[1, 1]
                    pooled_spec_logit <- summary_subset$coefficients[2, 1]
                    se_sens <- summary_subset$coefficients[1, 2]
                    se_spec <- summary_subset$coefficients[2, 2]
                    
                    pooled_sens_sub <- plogis(pooled_sens_logit)
                    pooled_spec_sub <- plogis(pooled_spec_logit)
                    
                    sens_ci_lower_sub <- plogis(pooled_sens_logit - 1.96 * se_sens)
                    sens_ci_upper_sub <- plogis(pooled_sens_logit + 1.96 * se_sens)
                    spec_ci_lower_sub <- plogis(pooled_spec_logit - 1.96 * se_spec)
                    spec_ci_upper_sub <- plogis(pooled_spec_logit + 1.96 * se_spec)
                ''')
                
                sens = ro.r('pooled_sens_sub')[0]
                spec = ro.r('pooled_spec_sub')[0]
                sens_ci = (ro.r('sens_ci_lower_sub')[0], ro.r('sens_ci_upper_sub')[0])
                spec_ci = (ro.r('spec_ci_lower_sub')[0], ro.r('spec_ci_upper_sub')[0])
                
                results[setting] = {
                    'sensitivity': sens,
                    'specificity': spec,
                    'sens_ci': sens_ci,
                    'spec_ci': spec_ci
                }
                
                print(f"Pooled Sensitivity: {sens:.3f} (95% CI: {sens_ci[0]:.3f} - {sens_ci[1]:.3f})")
                print(f"Pooled Specificity: {spec:.3f} (95% CI: {spec_ci[0]:.3f} - {spec_ci[1]:.3f})")
                
            except Exception as e:
                print(f"Warning: Could not perform meta-analysis for {setting}: {str(e)}")
                results[setting] = None
        
        return results
    
    def create_forest_plots(self, output_prefix='forest_plot'):
        """
        Create forest plots for sensitivity and specificity with subgroups
        """
        print("\n" + "="*70)
        print("CREATING FOREST PLOTS WITH SUBGROUPS")
        print("="*70)
        
        # Transfer data to R
        r_df = pandas2ri.py2rpy(self.df)
        ro.globalenv['meta_data'] = r_df
        
        for metric in ['Sensitivity', 'Specificity']:
            print(f"\nCreating forest plot for {metric}...")
            
            filename = f"{output_prefix}_{metric.lower()}.png"
            
            ro.r(f'''
                library(meta)
                library(metafor)
                
                # Prepare data for forest plot
                if ("{metric}" == "Sensitivity") {{
                    meta_data$prop <- meta_data$Sensitivity
                    meta_data$events <- meta_data$TP
                    meta_data$total <- meta_data$n_diseased
                    title_text <- "Forest Plot - Sensitivity by Setting"
                }} else {{
                    meta_data$prop <- meta_data$Specificity
                    meta_data$events <- meta_data$TN
                    meta_data$total <- meta_data$n_non_diseased
                    title_text <- "Forest Plot - Specificity by Setting"
                }}
                
                # Perform meta-analysis with subgroups using metaprop
                forest_meta <- metaprop(
                    event = events,
                    n = total,
                    studlab = Study_ID,
                    data = meta_data,
                    sm = "PLOGIT",
                    method = "Inverse",
                    method.tau = "DL",
                    byvar = Setting,
                    print.byvar = TRUE,
                    random = TRUE,
                    fixed = FALSE
                )
                
                # Create forest plot
                png("/mnt/user-data/outputs/{filename}", 
                    width=14, height=10, units="in", res=300)
                
                forest(forest_meta,
                       rightcols = c("effect", "ci"),
                       rightlabs = c("{metric}", "95% CI"),
                       leftcols = c("studlab", "event", "n"),
                       leftlabs = c("Study", "Events", "Total"),
                       xlab = "{metric}",
                       smlab = "",
                       weight.study = "random",
                       squaresize = 0.5,
                       col.square = "navy",
                       col.square.lines = "navy",
                       col.diamond = "maroon",
                       col.diamond.lines = "maroon",
                       pooled.totals = TRUE,
                       comb.fixed = FALSE,
                       comb.random = TRUE,
                       print.tau2 = TRUE,
                       print.I2 = TRUE,
                       digits = 3,
                       fontsize = 10,
                       spacing = 1.2,
                       addrow.subgroups = TRUE)
                
                title(title_text, line = 1, font.main = 2)
                
                dev.off()
            ''')
            
            print(f"  ✓ Saved: {filename}")
        
        print(f"\n✓ Forest plots created successfully!")
        return True
    
    def create_sroc_curve(self, output_file='sroc_curve.png'):
        """
        Create Summary ROC (SROC) curve
        """
        print("\n" + "="*70)
        print("CREATING SROC CURVE")
        print("="*70)
        
        # Transfer data to R
        r_df = pandas2ri.py2rpy(self.df)
        ro.globalenv['meta_data'] = r_df
        
        ro.r(f'''
            library(mada)
            
            # Fit bivariate model for SROC
            fit_sroc <- reitsma(meta_data, 
                               correction = 0.5,
                               correction.control = "all")
            
            # Create SROC plot with subgroups
            png("/mnt/user-data/outputs/{output_file}", 
                width=12, height=10, units="in", res=300)
            
            par(mar=c(5,5,4,2))
            
            # Get unique settings for colors
            settings <- unique(meta_data$Setting)
            colors <- c("blue", "red", "green", "purple")
            
            # Plot SROC curve
            plot(fit_sroc, 
                 sroclwd = 2,
                 main = "Summary ROC Curve with Confidence Region",
                 cex.main = 1.5,
                 xlim = c(0, 1),
                 ylim = c(0, 1),
                 xlab = "False Positive Rate (1 - Specificity)",
                 ylab = "True Positive Rate (Sensitivity)",
                 cex.lab = 1.2)
            
            # Add individual study points colored by setting
            for (i in 1:length(settings)) {{
                setting_data <- meta_data[meta_data$Setting == settings[i], ]
                fpr <- 1 - setting_data$Specificity
                tpr <- setting_data$Sensitivity
                
                points(fpr, tpr, 
                       pch = 19, 
                       col = colors[i],
                       cex = 1.5)
            }}
            
            # Add legend
            legend("bottomright", 
                   legend = settings,
                   col = colors[1:length(settings)],
                   pch = 19,
                   cex = 1.2,
                   title = "Setting",
                   bty = "n")
            
            # Add diagonal reference line
            abline(a = 0, b = 1, lty = 2, col = "gray", lwd = 1)
            
            # Add grid
            grid(col = "lightgray", lty = "dotted")
            
            dev.off()
        ''')
        
        print(f"✓ SROC curve created: {output_file}")
        return True
    
    def create_coupled_forest_plot(self, output_file='coupled_forest_plot.png'):
        """
        Create coupled forest plot showing sensitivity and specificity together
        """
        print("\n" + "="*70)
        print("CREATING COUPLED FOREST PLOT")
        print("="*70)
        
        r_df = pandas2ri.py2rpy(self.df)
        ro.globalenv['meta_data'] = r_df
        
        ro.r(f'''
            library(mada)
            
            # Create coupled forest plot
            png("/mnt/user-data/outputs/{output_file}", 
                width=16, height=10, units="in", res=300)
            
            par(mfrow=c(1,2), mar=c(5,4,4,2))
            
            # Sensitivity forest plot
            sens_meta <- metaprop(
                event = TP,
                n = n_diseased,
                studlab = Study_ID,
                data = meta_data,
                sm = "PLOGIT",
                method.tau = "DL",
                byvar = Setting,
                random = TRUE,
                fixed = FALSE
            )
            
            forest(sens_meta,
                   rightcols = c("effect", "ci"),
                   rightlabs = c("Sensitivity", "95% CI"),
                   leftcols = c("studlab"),
                   leftlabs = c("Study"),
                   xlab = "Sensitivity",
                   comb.fixed = FALSE,
                   comb.random = TRUE,
                   col.diamond = "maroon",
                   fontsize = 9)
            
            title("A) Sensitivity", adj = 0, font.main = 2)
            
            # Specificity forest plot
            spec_meta <- metaprop(
                event = TN,
                n = n_non_diseased,
                studlab = Study_ID,
                data = meta_data,
                sm = "PLOGIT",
                method.tau = "DL",
                byvar = Setting,
                random = TRUE,
                fixed = FALSE
            )
            
            forest(spec_meta,
                   rightcols = c("effect", "ci"),
                   rightlabs = c("Specificity", "95% CI"),
                   leftcols = c("studlab"),
                   leftlabs = c("Study"),
                   xlab = "Specificity",
                   comb.fixed = FALSE,
                   comb.random = TRUE,
                   col.diamond = "maroon",
                   fontsize = 9)
            
            title("B) Specificity", adj = 0, font.main = 2)
            
            dev.off()
        ''')
        
        print(f"✓ Coupled forest plot created: {output_file}")
        return True
    
    def create_summary_table(self, output_file='summary_table.csv'):
        """
        Create summary table of all studies with calculated metrics
        """
        print("\n" + "="*70)
        print("CREATING SUMMARY TABLE")
        print("="*70)
        
        summary_df = self.df[[
            'Study_ID', 'Setting', 'Device', 'Reference_Test',
            'Sensitivity', 'Specificity', 'AUC',
            'n_diseased', 'n_non_diseased',
            'TP', 'FN', 'FP', 'TN'
        ]].copy()
        
        # Calculate additional metrics
        summary_df['PPV'] = summary_df['TP'] / (summary_df['TP'] + summary_df['FP'])
        summary_df['NPV'] = summary_df['TN'] / (summary_df['TN'] + summary_df['FN'])
        summary_df['Accuracy'] = (summary_df['TP'] + summary_df['TN']) / (
            summary_df['TP'] + summary_df['TN'] + summary_df['FP'] + summary_df['FN']
        )
        
        # Save to outputs
        output_path = f'/mnt/user-data/outputs/{output_file}'
        summary_df.to_csv(output_path, index=False)
        
        print(f"✓ Summary table saved: {output_file}")
        print(f"\nTable preview:")
        print(summary_df.head())
        
        return summary_df
    
    def run_complete_analysis(self):
        """
        Run the complete meta-analysis workflow
        """
        print("\n" + "="*70)
        print("STARTING COMPLETE META-ANALYSIS WORKFLOW")
        print("="*70)
        
        # Step 1: Prepare data
        self.prepare_data()
        
        # Step 2: Overall meta-analysis
        overall_sens, overall_spec = self.perform_overall_meta_analysis()
        
        # Step 3: Subgroup meta-analysis
        subgroup_results = self.perform_subgroup_meta_analysis()
        
        # Step 4: Create forest plots
        self.create_forest_plots()
        
        # Step 5: Create SROC curve
        self.create_sroc_curve()
        
        # Step 6: Create coupled forest plot
        self.create_coupled_forest_plot()
        
        # Step 7: Create summary table
        summary_df = self.create_summary_table()
        
        print("\n" + "="*70)
        print("ANALYSIS COMPLETE!")
        print("="*70)
        print("\nGenerated files in /mnt/user-data/outputs/:")
        print("  - forest_plot_sensitivity.png")
        print("  - forest_plot_specificity.png")
        print("  - sroc_curve.png")
        print("  - coupled_forest_plot.png")
        print("  - summary_table.csv")
        print("\nAll outputs are ready for your living systematic review!")
        
        return {
            'overall_sensitivity': overall_sens,
            'overall_specificity': overall_spec,
            'subgroup_results': subgroup_results,
            'summary_table': summary_df
        }


# Main execution
if __name__ == "__main__":
    # Initialize meta-analysis with your CSV file
    meta = DiagnosticMetaAnalysis('niri_meta.csv')
    
    # Run complete analysis
    results = meta.run_complete_analysis()
    
    print("\n" + "="*70)
    print("FINAL SUMMARY")
    print("="*70)
    print(f"\nOverall Pooled Sensitivity: {results['overall_sensitivity']:.3f}")
    print(f"Overall Pooled Specificity: {results['overall_specificity']:.3f}")
    
    if results['subgroup_results']:
        print("\nSubgroup Results:")
        for setting, data in results['subgroup_results'].items():
            if data:
                print(f"\n  {setting}:")
                print(f"    Sensitivity: {data['sensitivity']:.3f}")
                print(f"    Specificity: {data['specificity']:.3f}")
    
    print("\n" + "="*70)
    print("Ready for your living systematic review!")
    print("="*70)
