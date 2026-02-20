# SUNFISH Part 2 Trial Simulation

This folder contains simulation scripts replicating the **SUNFISH Part 2** clinical trial for Bayesian subgroup analysis using the `bonsaiforest2` package.

## Trial Overview

**SUNFISH Part 2** was a randomized (2:1), double-blind, placebo-controlled trial evaluating risdiplam in patients with later-onset spinal muscular atrophy (SMA).

### Design Parameters
- **Sample Size**: N = 180 (120 Risdiplam : 60 Placebo)
- **Primary Endpoint**: Change from baseline in Motor Function Measure 32 (MFM32) total score at Month 12
- **Treatment Effect**: 1.55 points (Risdiplam: +1.36 vs Placebo: -0.19)
- **Standard Deviation**: 6.0 points (assumed for power calculations)

### Covariates

Based on Table 1, Supplementary Appendix Figure S3, and Forest Plot subgroup analysis:

1. **x_1 - Age Group** (Stratification factor)
   - 2-5 years: 31%
   - 6-11 years: 32%
   - 12-17 years: 25%
   - 18-25 years: 12%

2. **x_2 - SMA Type**
   - Type 2: 71%
   - Type 3: 29%

3. **x_3 - SMN2 Copy Number**
   - 2 copies: 3%
   - 3 copies: 87%
   - 4 copies: 10%

4. **x_4 - Disease Severity** (Based on baseline MFM32 quartiles)
   - Severe (≤Q1=37.5): ~25%
   - Moderate (>Q1 to Q3=54.17): ~50%
   - Mild (>Q3): ~25%

5. **x_5 - Region** (From Forest Plot subgroup analysis)
   - Europe: 67%
   - North America: 13%
   - China: 9.5%
   - Japan: 8.5%
   - Rest of World (RoW): 2%

6. **baseline_mfm** - Baseline MFM32 Score (continuous)
   - Mean: 46.11
   - SD: 11.46
   - Range: 0-96

## Simulation Scenarios

We simulate **6 scenarios** adapted from the original bonsaiforest2 framework to match the SUNFISH trial context:

### Scenario 1: Homogeneous Effect
Replicates the primary efficacy analysis result:
- Consistent **1.55 point improvement** across all subgroups
- Represents the overall population treatment effect reported in the trial

### Scenario 2: Region-Driven Heterogeneity
Replicates the dramatic regional differences observed in the Forest Plot subgroup analysis:
- **Europe**: **+2.40 points** (strong benefit, reference group)
- **North America**: **+2.44 points** (similar to Europe)
- **Japan**: **-3.58 points** (harm - qualitative interaction/sign reversal)
- **China**: **-1.85 points** (harm - qualitative interaction/sign reversal)
- **Rest of World**: **+3.30 points** (strong benefit)

This pattern demonstrates a **qualitative interaction** where treatment is beneficial in Western regions but harmful in Asian regions, representing the most extreme form of treatment heterogeneity.

### Scenario 3: Null Effect with Type-Based Crossover
Tests overall null treatment effect with SMA Type interaction:
- Overall mean treatment effect: **0 points**
- SMA Type 2 patients: **+2.0 points** (benefit)
- SMA Type 3 patients: **-2.0 points** (harm)
- Tests ability to detect crossover/qualitative interactions when overall effect is null

### Scenario 4: Mild Random Heterogeneity
Tests robustness to small, unpredictable heterogeneity:
- Base treatment effect: **1.55 points** (population average)
- Random noise (SD = 0.5) added to all subgroup interaction coefficients
- Represents realistic scenario where treatment effects vary slightly and unpredictably

### Scenario 5: Large Random Heterogeneity
Tests robustness to strong, unpredictable heterogeneity:
- Base treatment effect: **1.55 points** (population average)
- Large random noise (SD = 2.0) added to all subgroup interaction coefficients
- Stress test for shrinkage priors' ability to handle noisy, high-variance interactions

### Scenario 6: Age × SMA Type Interaction
Tests detection of complex 2-way interactions:
- Base effect: **2.5 points** (young Type 2 patients)
- Age effect: Progressive decline (-0.5 to -2.0 per age group)
- Type 3 effect: **-0.8 points**
- **Specific interaction**: Being both old (18-25y) AND Type 3 has additional penalty (**-1.5 points**)
- Tests ability to identify complex multiplicative interactions between categorical covariates

## File Structure

```
SUNFISH/
├── README.md                    # This file
├── functions.R                  # SUNFISH-specific simulation functions
├── Scenarios_generation.Rmd     # Script to generate datasets
├── Truth.Rmd                    # Script to compute true treatment effects
├── Results_analysis.Rmd         # Script to analyze simulation results
├── run_all.sh                   # Parallel job submission script
│
├── Scenarios/                   # Generated datasets (after running Scenarios_generation.Rmd)
│   ├── SUNFISH_Scenario_1.rds  # 1000 datasets - Homogeneous effect
│   ├── SUNFISH_Scenario_2.rds  # 1000 datasets - Age heterogeneity
│   ├── SUNFISH_Scenario_3.rds  # 1000 datasets - Null/Crossover
│   ├── SUNFISH_Scenario_4.rds  # 1000 datasets - Mild random heterogeneity
│   ├── SUNFISH_Scenario_5.rds  # 1000 datasets - Large random heterogeneity
│   └── SUNFISH_Scenario_6.rds  # 1000 datasets - Age×Type interaction
│
├── Results/                     # Analysis results (after running model scripts)
│
└── Model Scripts:
    ├── Horseshoe_low.R          # Global model with Horseshoe(scale_global=0.025)
    ├── Horseshoe_mid.R          # Global model with Horseshoe(scale_global=0.05)
    ├── Horseshoe_strong.R       # Global model with Horseshoe(scale_global=0.1)
    ├── R2D2_low.R               # Global model with R2D2(mean_R2=0.5, prec_R2=1, cons_D2=0.5)
    ├── R2D2_mid.R               # Global model with R2D2(mean_R2=0.5, prec_R2=2, cons_D2=1)
    ├── R2D2_strong.R            # Global model with R2D2(mean_R2=0.5, prec_R2=4, cons_D2=2)
    ├── Population.R             # Naive population (overall) analysis
    └── Subgroup.R               # Naive subgroup analysis
```

## Workflow

### 1. Generate Simulation Datasets

First, generate the 1000 replicate datasets for each scenario:

```r
# In R or RStudio
rmarkdown::render("Scenarios_generation.Rmd")
```

This creates:
- `Scenarios/SUNFISH_Scenario_1.rds` through `Scenarios/SUNFISH_Scenario_6.rds`
- Each file contains 1000 replicate datasets

### 2. Compute True Treatment Effects

Calculate the theoretical true treatment effects for all scenarios:

```r
# In R or RStudio
rmarkdown::render("Truth.Rmd")
```

This creates:
- `Scenarios/truth.RData` containing true overall and subgroup-specific treatment effects

### 3. Run Analysis Models

You can run individual model scripts:

```bash
# Run a single model
Rscript Horseshoe_low.R

# Or run all models in parallel (requires ACE platform)
./run_all.sh
```

Results are saved to the `Results/` folder.

### 4. Analyze Results

After all models complete, generate comprehensive performance analysis:

```r
# In R or RStudio
rmarkdown::render("Results_analysis.Rmd")
```

This produces:
- RMSE by scenario and estimator
- Standardized RMSE plots
- Bias and coverage analysis
- Subgroup-level performance comparisons
- Visual comparisons highlighting Scenario 2's regional heterogeneity

**Result files** in `Results/`:
- `SUNFISH_global_horseshoe_low.rds`
- `SUNFISH_global_horseshoe_mid.rds`
- `SUNFISH_global_horseshoe_strong.rds`
- `SUNFISH_global_r2d2_low.rds`
- `SUNFISH_global_r2d2_mid.rds`
- `SUNFISH_global_r2d2_strong.rds`
- `SUNFISH_population.rds`
- `SUNFISH_subgroup.rds`

## Key Differences from Standard Simulations

1. **Sample Size**: N=180 (vs N=1000 in standard simulations)
2. **Randomization**: 2:1 ratio (vs 1:1)
3. **Scenarios**: Only 2 scenarios (vs 6 in standard simulations)
   - Scenario 1: Homogeneous
   - Scenario 2: Age heterogeneity
4. **Covariates**: 4 specific SUNFISH covariates + baseline score (vs 10 generic)
5. **Effect Sizes**: Calibrated to SUNFISH trial results (1.55 overall, 3.14 to -0.65 by age)
6. **No OVAT Models**: OVAT (One-Variable-At-A-Time) not applicable with only 2 scenarios

## References

- Main Article: SUNFISH Part 2 clinical trial publication
- Supplementary Appendix: Figure S3 (subgroup analyses), Table 1 (baseline characteristics)

## Technical Notes

### Covariate Generation

The `simul_covariates_sunfish()` function generates correlated covariates matching the SUNFISH trial distribution. Key features:

#### Correlation Structure (from SAP and Trial Documentation)

**Base correlation**: All biological factors have moderate correlation (ρ = 0.2) reflecting shared disease biology.

**Specific correlations based on documented relationships**:

1. **Age × SMN2 Copy Number** (ρ = 0.3)
   - Type 2 SMA (more common in younger patients) often has fewer SMN2 copies
   - Stronger correlation reflects known biological relationship

2. **SMN2 Copy Number × Baseline MFM32** (ρ = 0.35)
   - More SMN2 copies generally associated with better motor function
   - Patients with 2 copies show decline; 3-4 copies show improvement (trial finding)

3. **Age × Baseline MFM32** (ρ = -0.15)
   - Older patients may have experienced more disease progression over time
   - Weak negative correlation reflects cumulative disease burden

4. **Disease Severity × All Variables** (deterministic)
   - Derived from baseline MFM32 quartiles (Q1=37.5, Q3=54.17)
   - Inherits all correlations from baseline MFM32

**Variable generation details**:
- Uses multivariate normal latents with specified correlation matrix
- Age, SMA Type, and SMN2 are categorized from normal quantiles
- Baseline MFM32 is continuous (Mean=46.11, SD=11.46, clamped to [0,96])
- Disease Severity derived from baseline MFM32 quartiles

### Model Formulation

The simulation uses the same formula structure as the main `bonsaiforest2` workflow:

```r
response_formula = "y ~ arm"
subgroup_model = ~ x_1 + x_2 + x_3 + x_4
```

Where:
- `y` = Change from baseline in MFM32
- `arm` = Treatment assignment (0=Placebo, 1=Risdiplam)
- Interactions between `arm` and subgroups capture heterogeneous treatment effects

### Prior Specifications

- **Horseshoe Priors**: Varying scale_global (0.025, 0.05, 0.1) for different shrinkage strengths
- **R2D2 Priors**: Varying precision parameters for alternative shrinkage approach
- All priors applied to predictive (treatment interaction) terms to identify heterogeneity

## Contact

For questions about the simulation setup, refer to the main `bonsaiforest2` package documentation or the development guide in `.github/copilot-instructions.md`.
