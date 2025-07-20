# Power Calculations for Parkinson’s and Metabolic Trait Analyses

This repository contains R code and a reproducible workflow to calculate power and required sample sizes for four core analyses in Parkinson’s disease (PD) and metabolic risk research.

## Analyses Covered

1. **GWAS (Case-Control Study)** – Estimate the total sample size needed to detect a genetic variant associated with PD.
2. **Phenotypic Clustering and Progression Modeling** – Determine the sample size for detecting small effects of biomarker-defined subtypes on disease progression.
3. **Mendelian Randomisation (MR)** – Estimate the sample size for detecting a causal effect of genetically predicted traits (e.g., insulin resistance) on PD.
4. **Propensity-Weighted Longitudinal Analysis** – Calculate the sample size needed to detect treatment effects (e.g., GLP-1RAs, ARBs) on PD outcomes in real-world datasets.

## Parameter Descriptions

| Parameter     | Description |
|---------------|-------------|
| `alpha`       | The significance level (typically 0.05), representing the probability of a Type I error (false positive). |
| `power`       | The statistical power (typically 0.80), representing the probability of correctly detecting a true effect (1 − Type II error). |
| `OR`          | The odds ratio representing the effect size in GWAS. OR > 1 indicates increased risk; OR < 1 indicates protection. |
| `MAF`         | Minor Allele Frequency, the frequency of the less common allele in the population; used in GWAS calculations. |
| `h`           | Cohen’s h, a standardized effect size used for comparing two proportions (e.g., allele frequency differences). |
| `f2_lmm`      | Cohen’s f², a standardized effect size for multiple regression or linear mixed models. f² = R² / (1 − R²). |
| `u_lmm`       | The number of predictor variables or degrees of freedom in the numerator (e.g., number of biomarker axes). |
| `r2_mr`       | The proportion of variance in the exposure explained by the genetic instrument in Mendelian randomisation. |
| `d`           | Cohen’s d, a standardized mean difference used for two-sample comparisons (e.g., treatment vs. control). |

## Usage

Run the R script `power_calculations.R` to print sample size estimates for each analysis to the console.

```r
install.packages("pwr")
source("power_calculations.R")
```
