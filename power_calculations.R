# Power Calculations for Parkinson’s Disease and Metabolic Trait Analyses
# Author: Your Name
# Date: YYYY-MM-DD

# Load required packages
library(pwr)

# Set common parameters
alpha <- 0.05
power <- 0.8

# ==========================================
# 1. GWAS (Case-Control Study)
# Objective: Detect a common variant (MAF = 0.2) associated with PD (OR = 1.3)
# Required sample size: 2,644 (1,322 PD cases, 1,322 controls)
# Recommended cohorts:
#   - UK Biobank (UKBB): >4,000 PD, >35,000 T2D, with genotype and metabolic traits
#   - All of Us: >3,000 PD, >30,000 T2D, ancestry diversity, genotyped
#   - Our Future Health: >2,000 PD, 35,000 T2D
#   - PPMI: ~2,000 PD, subset with T2D (~100–200)
#   - OPDC: ~1,000 PD, limited T2D (~<50 cases)
#   - Tracking Parkinson’s: ~1,800 PD, estimated ~90–100 T2D cases
# ==========================================

OR <- 1.3
MAF <- 0.2

p1 <- MAF
odds1 <- p1 / (1 - p1)
odds2 <- odds1 * OR
p2 <- odds2 / (1 + odds2)
h <- 2 * asin(sqrt(p1)) - 2 * asin(sqrt(p2))

gwas_power <- pwr.2p.test(h = h, sig.level = alpha, power = power, alternative = "two.sided")
n_gwas_total <- ceiling(gwas_power$n * 2)
gwas_per_group <- ceiling(gwas_power$n)
cat("GWAS Total N:", n_gwas_total, "| Per Group:", gwas_per_group, "\n\n")

# ==========================================
# 2. Phenotypic Clustering + Disease Progression (LMM / Cox)
# Objective: Test whether biomarker-defined PD subtypes predict clinical progression
# Required sample size: 647 (to detect small effect size, f² = 0.02, 5 predictors)
# Recommended cohorts:
#   - PPMI: ~2,000 PD, SomaScan, subset T2D (~100–200)
#   - UKBB (subset with Olink): >50,000 proteomes, deep metabolic traits
#   - Tracking Parkinson’s: 1,800 PD, SomaScan; ~90–100 T2D
#   - OPDC: ~1,000 PD with Olink, <50 T2D
#   - GNPC: >8,000 PD with harmonized SomaScan proteomics across multiple cohorts
# ==========================================

f2_lmm <- 0.02
u_lmm <- 5

lmm_power <- pwr.f2.test(u = u_lmm, f2 = f2_lmm, sig.level = alpha, power = power)
n_lmm_total <- ceiling(lmm_power$u + lmm_power$v + 1)
cat("LMM/Cox Total N:", n_lmm_total, "\n\n")

# ==========================================
# 3. Propensity-Weighted Longitudinal Model (Drug Repurposing)
# Objective: Evaluate treatment effect of GLP-1R agonists or ARBs in PD/T2D individuals
# Required sample size: 199 (≈100 per group; Cohen's d = 0.4 ~ HR = 0.7)
# Recommended cohorts:
#   - CPRD: >250,000 PD, >500,000 T2D; detailed RAS/GLP-1RA drug exposure
#   - UKBB: Exploratory use; fewer exposed individuals but linked outcomes
# ==========================================

d_propensity <- 0.4
prop_power <- pwr.t.test(d = d_propensity, power = power, sig.level = alpha, 
                         type = "two.sample", alternative = "two.sided")
n_propensity_total <- ceiling(prop_power$n * 2)
n_propensity_per_group <- ceiling(prop_power$n)
cat("Propensity-Weighted Model Total N:", n_propensity_total, 
    "| Per Group:", n_propensity_per_group, "\n")
