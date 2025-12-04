# Count Data Simulation Parameters - Analysis & Recommendation

## Problem

The count data simulations are producing "weird results" because they have **much less variance** than the TTE simulations, leading to very different statistical power:

- **TTE**: N=1000, events=247, Power ≈ 90%, SE(log HR) = 0.1273
- **Count (current)**: N=650, φ=2, Power ≈ 100%, SE(log RR) = 0.077

The count data has only **37%** of the variance of TTE data, making it much easier to detect effects and not comparable for simulation studies.

## Root Cause

For negative binomial distribution:
- Variance = μ + μ²/φ
- Current overdispersion φ=2 is too high (less variance)
- With smaller sample size (650 vs 1000), count still has higher power

## Recommended Solution

**Change count simulation parameters in `functions.R`:**

```r
"count" = list(
  endpoint = "count",
  N = 650,              # Keep current N
  intercept = 1.0,      
  base_effect = log(0.66),
  overdispersion = 0.35  # CHANGE from 2.0 to 0.35
)
```

### Justification

This gives:
- **Power ≈ 83%** (closer to TTE's 90%, and above target 80%)
- **SE(log RR) ≈ 0.143** (comparable to TTE's 0.127)
- More realistic overdispersion (typical count data has φ < 1)

### Alternative Options Tested

| N    | φ    | Power  | SE     | Notes                              |
|------|------|--------|--------|------------------------------------|
| 1000 | 0.50 | 98.7%  | 0.0992 | Match TTE N, still too much power  |
| 650  | 0.35 | 82.8%  | 0.1429 | **RECOMMENDED - balanced**         |
| 400  | 1.00 | 93.0%  | 0.1209 | Moderate, slightly high power      |
| 250  | 2.00 | 91.7%  | 0.1241 | Very small N, less stable          |

## Implementation

Update line 231 in `/home/pedreram/bonsaiforest2/simulations/functions.R`:

```r
overdispersion = 0.35  # Changed from 2.0
```

Then regenerate all count scenarios and rerun simulations.

## Technical Details

### Negative Binomial Variance Structure
- Mean: μ = exp(η), where η is linear predictor
- Variance: Var(Y) = μ + μ²/φ
- Lower φ → higher variance → more realistic count data
- φ ≈ 0.35 gives coefficient of variation ≈ 1.8 (high but realistic for overdispersed counts)

### Power Calculation for Count Data
- SE(log RR) ≈ √[Var(Y_C)/(n_C·μ_C²) + Var(Y_T)/(n_T·μ_T²)]
- With RR = 0.66, two-sided α = 0.05
- Target power 80% requires SE ≈ 0.14-0.15

## Verification: All Endpoints

After fixing Count, here's the final comparison across all 4 endpoints:

| Endpoint   | N    | SE     | Power | Variance Ratio | Status |
|------------|------|--------|-------|----------------|--------|
| TTE        | 1000 | 0.1273 | 90.4% | 1.00 (ref)     | ✓      |
| Binary     | 760  | 0.1490 | 79.7% | 1.37           | ✓      |
| Count      | 650  | 0.1429 | 82.8% | 1.26           | ✓      |
| Continuous | 350  | 0.1069 | 80.1% | 0.71           | ✓      |

**Assessment**: All endpoints now have comparable power (80-90%) and variance ratios within acceptable range (0.7-1.4x). The variance differences reflect the natural statistical properties of each endpoint type and are appropriate for simulation comparisons.

**Key Achievement**: Count variance ratio improved from 0.37x → 1.26x, making it properly comparable to TTE.

---
**Date**: November 27, 2025
**Issue**: Count simulations not comparable to TTE due to variance mismatch
**Status**: ✅ RESOLVED - All endpoints verified and balanced

---
**Date**: November 28, 2025
**Additional Issue**: Scenarios 4 & 5 had different heterogeneity across endpoints
**Fix Applied**: Scaled Binary/Count/Continuous coefficients by 0.85 in `functions.R`
**Status**: ✅ RESOLVED - All endpoints now have identical heterogeneity (SD=0.1248 for Scenario 4, SD=0.2496 for Scenario 5)
