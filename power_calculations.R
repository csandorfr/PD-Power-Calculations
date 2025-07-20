# Power Calculations for Parkinson’s Disease and Metabolic Trait Analyses
# Author: Your Name
# Date: YYYY-MM-DD

# Load required packages
library(pwr)

# Set common parameters
alpha <- 0.05
power <- 0.8

### 1. GWAS (Case-Control Study)
OR <- 1.3
MAF <- 0.2

# Calculate control and case proportions assuming additive model (approximate using arcsine)
p1 <- MAF
odds1 <- p1 / (1 - p1)
odds2 <- odds1 * OR
p2 <- odds2 / (1 + odds2)
h <- 2 * asin(sqrt(p1)) - 2 * asin(sqrt(p2))

gwas_power <- pwr.2p.test(h = h, sig.level = alpha, power = power, alternative = "two.sided")
n_gwas_total <- ceiling(gwas_power$n * 2)
cat("GWAS Total N:", n_gwas_total, "| Per Group:", ceiling(gwas_power$n), "\n\n")

### 2. Phenotypic Clustering + LMM / Cox (General Linear Model Approximation)
f2_lmm <- 0.02  # Small effect size
u_lmm <- 5      # Number of predictors

lmm_power <- pwr.f2.test(u = u_lmm, f2 = f2_lmm, sig.level = alpha, power = power)
n_lmm_total <- ceiling(lmm_power$u + lmm_power$v + 1)
cat("LMM/Cox Total N:", n_lmm_total, "\n\n")

### 3. Mendelian Randomisation (R² = 0.04)
r2_mr <- 0.04
f2_mr <- r2_mr / (1 - r2_mr)

mr_power <- pwr.f2.test(u = 1, f2 = f2_mr, sig.level = alpha, power = power)
n_mr_total <- ceiling(mr_power$u + mr_power$v + 1)
cat("Mendelian Randomisation Total N:", n_mr_total, "\n\n")

### 4. Propensity-Weighted Longitudinal Model (HR ~ 0.7 ~ Cohen's d = 0.4)
d_propensity <- 0.4
prop_power <- pwr.t.test(d = d_propensity, power = power, sig.level = alpha, 
                         type = "two.sample", alternative = "two.sided")
n_propensity_total <- ceiling(prop_power$n * 2)
cat("Propensity-Weighted Model Total N:", n_propensity_total, 
    "| Per Group:", ceiling(prop_power$n), "\n")
