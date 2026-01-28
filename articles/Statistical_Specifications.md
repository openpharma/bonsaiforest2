# Statistical Specifications

## 1 Scope of this document

This document describes the statistical methods implemented in the
`bonsaiforest2` R package for Bayesian subgroup analysis of continuous,
binary, count, and time-to-event outcomes in randomized clinical trials
(RCTs).

The package unable the implementation of both a **global modeling
approach** ([Wolbers et al. 2025](#ref-wolbers2025using)), that
estimates all main (prognostic) and interaction (predictive) effects in
a single model, and a **One-variable-at a time** ([Wang et al.
2024](#ref-wang2024bayesian)), that estimates the interaction effects
using a different model for each subgrouping variable.

The core features described in this document are:

1.  **Flexible Model Specification:** A unified framework to distinguish
    between:
    - **Prognostic** (main) vs. **Predictive** (interaction) effects.
    - **Shrunk** vs. **Unshrunk** coefficients. Users can specify
      unshrunk terms (with weakly informative priors), and optionally
      add shrunk prognostic effects and/or shrunk predictive effects
      (with strong regularization priors) as needed for their analysis.
    - **Customizable Priors:** Users can specify custom priors for
      individual coefficients or coefficient groups when domain
      expertise is available. If no custom priors are specified, the
      package applies automatic, well-calibrated default priors based on
      established best practices.
2.  **Advanced Shrinkage Priors:** Supports state-of-the-art shrinkage
    priors, including the **Regularized Horseshoe** ([Piironen and
    Vehtari 2017](#ref-piironen2017sparsity); [Wolbers et al.
    2025](#ref-wolbers2025using)) and **R2D2** ([Zhang et al.
    2022](#ref-zhang2022bayesian)), to provide robust control of
    false-positive findings while identifying potential heterogeneity.
3.  **Marginal Effect Estimation:** Use of **standardization
    (G-computation)** to derive interpretable *marginal* subgroup
    treatment effects, which average over the covariate distributions
    within each subgroup ([Wolbers et al. 2025](#ref-wolbers2025using)).
    This is implemented through:
    - [`estimate_subgroup_effects()`](https://openpharma.github.io/bonsaiforest2/reference/estimate_subgroup_effects.md)
      for the computational G-computation workflow.
    - [`summary_subgroup_effects()`](https://openpharma.github.io/bonsaiforest2/reference/summary_subgroup_effects.md)
      for user-friendly effect summaries.

This document is structured as follows: [Section 2](#sec-intro) provides
a conceptual introduction to subgroup analysis challenges and the
Bayesian shrinkage solution. [Section 3](#sec-methodology) details the
statistical methodology, including notation, the global and the
one-variable at a time models, endpoint-specific likelihoods, prior
specifications, and the standardization procedure. [Section
4](#sec-functions) maps these statistical methods to the core functions
in the `bonsaiforest2` package.

## 2 Introduction to Subgroup Analysis and Estimation

### 2.1 The Challenge of Subgroup Analysis

Exploratory subgroup analyses are a standard part of RCT reporting,
typically visualized in forest plots to assess the consistency of
treatment effects. However, their interpretation is notoriously
difficult due to several statistical challenges:

- **Low Power:** Subgroups have smaller sample sizes, leading to low
  precision (wide confidence intervals) and an inability to detect true,
  modest differences in effect.
- **Multiplicity:** Examining many subgroups (e.g., 20 or more) inflates
  the risk of false-positive findings. Investigators may focus on the
  most extreme results, which are often spurious.
- **Overfitting:** Standard subgroup-specific estimates (fitting a model
  only to data within that subgroup) are unbiased but have high
  variance, leading to an overestimation of heterogeneity.

This has led to the common recommendation that the *overall* treatment
effect is often a more reliable estimate for a subgroup than the
estimate derived from the subgroup’s data alone.

### 2.2 Prognostic vs. Predictive Effects

To build a meaningful model, it is crucial to distinguish how a baseline
variable relates to the outcome:

- **Prognostic Effect**: A variable is **prognostic** if it is
  associated with the clinical outcome, regardless of the treatment
  received. It describes the natural course of the disease. In a model,
  this is a **main effect**.
- **Predictive Effect**: A variable is **predictive** if its value
  modifies the magnitude or direction of the treatment’s effect. It
  identifies heterogeneity. In a model, this is a
  **treatment-by-variable interaction**.

`bonsaiforest2` is built to model both types of effects explicitly and
apply different modeling strategies (e.g., shrinkage) to each.

### 2.3 The Role of Bayesian Shrinkage

Bayesian shrinkage methods directly address the challenges of subgroup
analysis by providing a principled compromise between the two extremes
of “no pooling” (subgroup-specific estimates) and “full pooling”
(overall population estimate).

This compromise, often called **partial pooling**, is achieved using
hierarchical models. Subgroup effects are assumed to be *exchangeable*
(drawn from a common distribution). This assumption allows subgroups to
“borrow strength” from each other:

- Subgroups with **small sample sizes** and extreme results are “shrunk”
  heavily toward the overall mean, reducing their variance.
- Subgroups with **large sample sizes** and strong evidence of a real
  effect are shrunk less, allowing their data to dominate.

This results in estimates with a better **bias-variance trade-off**: we
accept a small amount of bias (by pulling real effects slightly toward
the mean) in exchange for a large reduction in variance, leading to a
lower overall error and better control of spurious findings.

### 2.4 Conditional vs. Marginal Estimands

When covariates are included in a non-linear model (like a logistic or
Cox model), the resulting coefficients represent **conditional**
effects. For example, the treatment effect is the effect for a patient
with a *specific set* of covariate values (e.g., the reference values).

This is problematic for forest plots, where we want to know the
**marginal** effect: “What is the average treatment effect for *all
patients* in the ‘Female’ subgroup, averaging over their different ages,
regions, etc.?”.

Because different subgroups have different covariate distributions
(e.g., the `age >= 65` subgroup will have a different health status
profile than the `age < 65` group), these conditional coefficients are
not directly comparable.

`bonsaiforest2` solves this by using **Standardization (G-computation)**
([**wolbers2024shrinkage?**](#ref-wolbers2024shrinkage)). This procedure
correctly averages over the specific covariate distribution of each
subgroup to estimate a valid marginal treatment effect, ensuring all
subgroups are compared on the same interpretable scale. (See [Section
3.8](#sec-standardization) for a detailed explanation of how this
standardization is performed in the package.)

## 3 Statistical Methodology

### 3.1 Setting and Notation

We establish the following notation for the global model:

- \\i = 1, \dots, N\\: Index for individual patients.
- \\y_i\\: The outcome for patient \\i\\. We will consider that this
  variable can only be a binary, count, continuous or time-to-event
  endpoint.
- \\s_i\\: The binary treatment indicator (\\s_i=1\\ for treatment,
  \\s_i=0\\ for control).
- \\l = 1, \dots, L\\: Index for an overall subgroup level (e.g.,
  “Male”, “Female”, “EU”, “USA”).
- \\x\_{il}\\: Indicator (0/1) for patient \\i\\ belonging to subgroup
  level \\l\\ (used for prognostic effects).
- \\z\_{il} = s_i \cdot x\_{il}\\: The interaction between treatment and
  subgroup level \\l\\ (used for predictive effects).
- \\u\_{iv}\\: Value of the \\v\\-th additional non-subgroup covariate
  for patient \\i\\. This additional covariates don’t need to be binary
  or categorical.
- \\\boldsymbol{\mu}\\: The \\N \times 1\\ vector of expected outcomes.
- \\g(\cdot)\\: A link function. In [Section x](#sec-endpoint) we define
  the used link function for each type of endpoint.
- \\\boldsymbol{\eta}\\: The \\N \times 1\\ vector of the linear
  predictor.

### 3.2 Design Matrix Considerations

A critical implementation detail is the construction of the design
matrix for the subgroup-by-treatment interactions.

The package treats prognostic effects and predictive effects
differently, as they serve different purposes in the model:

- **For prognostic effects:** The goal is interpretability and avoiding
  collinearity with the intercept. Here, standard treatment (dummy)
  coding is used. For each of the J subgrouping variables, one level is
  treated as reference, resulting in L−J total columns for the subgroup
  main effects.

- **For predictive effects:** The package creates explicit dummy
  variables for each treatment-by-subgroup combination. For a
  subgrouping variable with levels {A, B, C}, the package creates
  interaction dummies `trt_subgroupA`, `trt_subgroupB`, `trt_subgroupC`,
  where each dummy equals 1 when the patient is both in the treatment
  group AND in that specific subgroup level. This approach ensures:

  - All subgroup levels are treated symmetrically under the
    exchangeability assumption
  - Each subgroup has an interpretable parameter directly representing
    its treatment effect modification
  - Predictions work correctly through G-computation (see Section 3.8)

### 3.3 The Principle of Model Hierarchy

A fundamental concept in regression modeling is the **principle of model
hierarchy** (or marginality). This principle states that if a model
includes an interaction term (e.g., `A:B`), it should also include the
corresponding main effects (e.g., `A` and `B`).

In the context of `bonsaiforest2`, this means that any subgrouping
variable included as **predictive** (i.e., in an interaction with
treatment, like `trt:subgroup`) *must* also be included as
**prognostic** (i.e., as a main effect, `subgroup`).

Including the main effect ensures that the interaction term purely
captures the *modification* of the treatment effect, rather than a mix
of the main prognostic effect and the interaction, leading to a more
stable and interpretable model.

The `bonsaiforest2` package enforces this principle automatically. The
`prepare_formula_model` function checks if a main effect is specified as
prognostic for every variable included in a predictive term. If it is
not, the function automatically adds the main effect to the list of
unshrunk prognostic terms to ensure the model hierarchy is respected.

### 3.4 The Global Model

The package enables the implementation of a global model ([Wolbers et
al. 2025](#ref-wolbers2025using)) where all effects are estimated in a
single, comprehensive regression formula.

#### 3.4.1 Full Model Equation

***Matrix form:***

The model for the linear predictor, \\\boldsymbol{\eta}\\, for all \\N\\
patients is given by:

\\g(\boldsymbol{\mu}) = \boldsymbol{\eta} =
\underbrace{\mathbf{X_p}\boldsymbol{\beta\_{1,p}}}\_{\text{Penalized
prognostic effects}}
+\underbrace{\mathbf{X_n}\boldsymbol{\beta\_{1,n}}}\_{\text{Non-penalized
prognostic effects}}
+\underbrace{\mathbf{Z_p}\boldsymbol{\beta\_{2,p}}}\_{\text{Penalized
predictive effects}} +\\ \\
\underbrace{\mathbf{Z_n}\boldsymbol{\beta\_{2,n}}}\_{\text{Non-penalized
predictive effects}}+
\underbrace{\mathbf{U}\boldsymbol{\beta\_{3}}}\_{\text{Extra prognostic
effects}} \\

The Bayesian hierarchical structure is applied as follows:

- **Penalized Coefficients**: \\\boldsymbol{\beta\_{1,p}}\\ and
  \\\boldsymbol{\beta\_{2,p}}\\ are given shrinkage priors.

- **Non-Penalized Coefficients**: \\\boldsymbol{\beta\_{1,n}}\\ and
  \\\boldsymbol{\beta\_{2,n}}\\ are given standard, weakly informative
  priors.

- **Extra prognostic effects Coefficients**: \\\boldsymbol{\beta_3}\\
  are given standard, weakly informative priors.

**Components definition:**

- \\\mathbf{X}= \[\mathbf{X_p} \\\\ \mathbf{X}\_n\]\\ is the \\N \times
  M\\ design matrix for the **prognostic effects**. It typically
  includes:

  - An intercept column. In [Section x](#sec-intercept) we discuss the
    inclusion and the prior we give to the intercept in detail

  - A column for the main treatment effect.

  - Columns for the prognostic effects of each subgroup level (using
    dummy coding, so one reference level per factor is dropped). Each
    prognostic effect will be contained in \\\mathbf{X}\_n\\ or
    \\\mathbf{X}\_p\\ depending if we want to shrink or not the effect.

    Example:

  ``` r
  # 1. Create your data as a data frame
  design_matrix_df <- data.frame(
    Patient    = c(1, 2, 3, 4),
    Intercept  = c(1, 1, 1, 1),
    trt        = c(0, 1, 0, 1),
    sexF       = c(0, 1, 1, 0),
    regionUS   = c(0, 0, 1, 1),
    regionAsia = c(0, 1, 0, 0)
  )

  # 2. Use knitr::kable() to format and print the table
  knitr::kable(
    design_matrix_df,
    caption = "An example design matrix for a global subgroup model.",
    align = c('l', 'c', 'c', 'c', 'c') # Align columns (l=left, c=center)
  )
  ```

  | Patient | Intercept | trt | sexF | regionUS | regionAsia |
  |:--------|:---------:|:---:|:----:|:--------:|:-----------|
  | 1       |     1     |  0  |  0   |    0     | 0          |
  | 2       |     1     |  1  |  1   |    0     | 1          |
  | 3       |     1     |  0  |  1   |    1     | 0          |
  | 4       |     1     |  1  |  0   |    1     | 0          |

  Table 3.1: An example design matrix for a global subgroup model.

  *Note:* M=1+1+L-J= 1 degree for the intercept+ 1 for the treatment
  effect + L (Number of subgroups) - J(Number of variables) because we
  get the reference levels out.

- \\\boldsymbol{\beta\_{1,p}}\\ is the \\P \times 1\\ vector of
  coefficients for the prognostic effects that we want to shrink. We
  will give an shrinkage prior to this coefficients.

- \\\boldsymbol{\beta\_{1,n}}\\ is the \\(M-P) \times 1\\ vector of
  coefficients for the prognostic effects that we don’t want to shrink.

- \\\mathbf{Z}= \[\mathbf{Z_p} \\\\ \mathbf{Z}\_n\]\\ is the \\N \times
  L\\ design matrix for the **predictive effects**. Its columns
  represent explicit dummy variables for the interaction between
  treatment and each of the \\L\\ subgroup levels. Each column
  `trt_subgroupLEVEL` equals 1 when patient \\i\\ is both on treatment
  (\\s_i = 1\\) AND belongs to that specific subgroup level, and 0
  otherwise. This ensures all subgroups are treated symmetrically.

``` r
# Manually creating the interaction matrix for clarity
# Each column represents trt_subgroupLEVEL (1 if treated AND in that level, 0 otherwise)
interaction_matrix <- data.frame(
  Patient           = 1:6,
  `trt_sexM`        = c(1, 0, 1, 1, 0, 0),
  `trt_sexF`        = c(0, 1, 0, 0, 1, 1),
  `trt_regionEU`    = c(0, 1, 0, 0, 1, 0),
  `trt_regionUS`    = c(1, 0, 0, 1, 0, 0),
  `trt_regionAsia`  = c(0, 0, 1, 0, 0, 1)
)

# Using knitr::kable for a nice output
knitr::kable(
  interaction_matrix,
  align = 'c',
  col.names = c("Patient", "trt_sexM", "trt_sexF", "trt_regionEU", "trt_regionUS", "trt_regionAsia"),
  caption = "Interaction design matrix using explicit dummy variables"
)
```

| Patient | trt_sexM | trt_sexF | trt_regionEU | trt_regionUS | trt_regionAsia |
|:-------:|:--------:|:--------:|:------------:|:------------:|:--------------:|
|    1    |    1     |    0     |      0       |      1       |       0        |
|    2    |    0     |    1     |      1       |      0       |       0        |
|    3    |    1     |    0     |      0       |      0       |       1        |
|    4    |    1     |    0     |      0       |      1       |       0        |
|    5    |    0     |    1     |      1       |      0       |       0        |
|    6    |    0     |    1     |      0       |      0       |       1        |

Table 3.2: Interaction design matrix using explicit dummy variables

- \\\boldsymbol{\beta\_{2,p}}\\ is the \\R \times 1\\ vector of
  coefficients for the predictive effects that we want to shrink. We
  will give an shrinkage prior to this coefficients.
- \\\boldsymbol{\beta\_{2,n}}\\ is the \\(L-R) \times 1\\ vector of
  coefficients for the predictive effects that we don’t want to shrink.
- \\U\\ is the \\N \times V\\ design matrix for the extra variables we
  want to take into account in the fitting of the model. This will also
  be **unpenalized.**
- \\\boldsymbol{\beta_3}\\ is the \\V \times 1\\ vector of extra
  prognostic effects coefficients

|                       | **Shrinkage**                             | **Not Shrinkage**                         |
|:----------------------|:------------------------------------------|:------------------------------------------|
| **Prognostic Effect** | \\\mathbf{X_p}\boldsymbol{\beta\_{1,p}}\\ | \\\mathbf{X_n}\boldsymbol{\beta\_{1,n}}\\ |
| **Predictive Effect** | \\\mathbf{Z_p}\boldsymbol{\beta\_{2,p}}\\ | \\\mathbf{Z_n}\boldsymbol{\beta\_{2,n}}\\ |

***Linear form:***

To make it more concrete, the linear predictor for a single patient
\\i\\ can be written as:

\\ \eta_i = \underbrace{\beta_0 + \beta\_{\text{treat}}s_i +
\underbrace{\sum\_{l=1}^{P} \beta^l\_{1,p}x^{il}\_n}\_\text{penalized}+
\underbrace{\sum\_{l=P+1}^{M} \beta^l\_{1,n}x^{il}\_n}\_\text{Not
penalized}+\underbrace{\sum\_{v=1}^{V} \beta\_{3}^vu^{iv}}\_{\text{Extra
}}}\_{\text{Prognostic Effects}} +
\underbrace{\underbrace{\sum\_{l=1}^{L-R} \beta\_{2,p}^l
z_p^{il}}\_{\text{penalized }} + \underbrace{\sum\_{l=L-R+1}^{L}
\beta\_{2,n}^l z_n^{il}}\_{\text{Not penalized}} }\_{\text{Predictive
effects}} \\

**Note**: See that even if \\z^{il}=x^{il}\cdot s^{i}\\, as the
prognostic effect or the predictive effect of the same subgrouping
variable may be penalized or not, we separate it in \\x\\ and \\z\\ in
the formula.

#### 3.4.2 Example

A randomized, double-blind trial comparing a new drug, **DrugX**,
against a **placebo** in patients with hypertension.

- **Primary Endpoint (**\\Y\_{i}\\): Change from baseline in systolic
  blood pressure (SBP) in mmHg at 6 months.
- **Subgroups of Interest**:
  - **Age Group**: `<65` vs. `≥65` years.
  - **Sex**: `Male` vs. `Female`.
  - **Kidney Function**: eGFR `<60` (impaired) vs. `≥60` (normal).
- **Covariates**: `BaselineSBPᵢ`, `BMIᵢ`, `Smokerᵢ` (binary).

**Global Model with Shrinkage**

The outcome is assumed to be normally distributed: \\Y\_{i} \sim
N(\mu\_{i}, \sigma^2)\\.

The linear predictor \\\mu\_{i}\\ is modeled as:

\\ \mu\_{i} = \underbrace{\beta_0}\_{\text{Intercept}} +
\underbrace{\beta\_{\text{treat}} \cdot
\text{Treatment}\_i}\_{\text{Main Trt Effect}} +
\underbrace{\sum\_{k=1}^{K} \alpha_k \cdot
\text{Subgroup}\_{ik}}\_{\text{Prognostic Subgroup Effects}} +
\underbrace{\sum\_{k=1}^{K} \gamma_k \cdot (\text{Subgroup}\_{ik} \cdot
\text{Treatment}\_i)}\_{\text{Shrunken Predictive Effects}} +
\underbrace{\sum\_{j=1}^{J} \delta_j \cdot
\text{Covariate}\_{ij}}\_{\text{Extra Prognostic effects}} \\

\\ =\underbrace{\beta_0}\_{\text{Intercept}} +
\underbrace{\beta\_{\text{treat}} \cdot
\text{Treatment}\_i}\_{\text{Main Trt Effect}} + \underbrace{\alpha_1
\cdot \text{Age Group}\_i+\alpha_2 \cdot \text{Sex}\_i+\alpha_3 \cdot
\text{eGFR}\_i}\_{\text{Prognostic Subgroup Effects}}+ \\

\\ + \underbrace{\gamma_1 \cdot (\text{Age Group}\_i\cdot
\text{Treatment}\_i)+\gamma_2 \cdot (\text{Sex}\_i\cdot
\text{Treatment}\_i)+\gamma_3 \cdot (\text{eGFR}\_i\cdot
\text{Treatment}\_i)}\_{\text{Predictive Effects}} +
\underbrace{\delta_1 \cdot \text{BaselineSBP}\_i+\delta_2 \cdot
\text{BMI}\_i+\delta_3 \cdot \text{Smoker}\_i}\_{\text{Extra Prognostic
effects}} \\

### 3.5 One-Variable-at-a-Time Model

The **one-variable-at-a-time model** ([Wang et al.
2024](#ref-wang2024bayesian)) is a Bayesian hierarchical model (BHM)
that improves upon standard subgroup analysis by borrowing information
across subgroup levels to produce more stable and reliable estimates.
Instead of analyzing each subgroup level in isolation, it analyzes all
levels *within a single subgrouping variable* (e.g., all age groups)
together in one model.

For example, for the variable *Sex* (with levels Male and Female), you
would fit one hierarchical model. In this single model, the treatment
effects for Male (\\\beta\_{2,\text{sex,male}}\\) and Female
(\\\beta\_{2,\text{sex,female}}\\) are treated as exchangeable and are
linked by a common prior distribution. Then, a completely separate
hierarchical model would be fit for the variable *Region*.

#### 3.5.1 Full Model Equation

For a given subgrouping variable \\j\\, we can specify a flexible linear
predictor that allows for some covariate effects to be constant across
the levels of \\j\\, while others are allowed to vary.

Let’s partition the patient-level covariates into two sets:

- \\\mathbf{u}\_i\\: A vector of covariates assumed to have a **constant
  prognostic effect** across all levels of the subgrouping variable
  \\j\\.

- \\\mathbf{v}\_i\\: A vector of covariates for which there is a strong
  *a priori* reason to believe their prognostic effect **varies** across
  the levels \\k\\ of the subgrouping variable \\j\\.

The linear predictor for patient \\i\\ in level \\k\\ of variable \\j\\
is then: \\g(\mu\_{i,j,k}) = \beta\_{1,j,k} + \beta\_{2,j,k}z_i +
\boldsymbol{\beta}\_{3,j}\mathbf{u}\_i +
\boldsymbol{\beta}\_{4,j,k}\mathbf{v}\_i\\

- \\\beta\_{2,j,k}\\: The predictive effect for level \\k\\ of variable
  \\j\\. This is the primary parameter of interest for shrinkage.
- \\\boldsymbol{\beta}\_{3,j}\mathbf{u}\_i\\: The effect of covariates
  assumed to be **constant** across levels. The coefficient vector
  \\\boldsymbol{\beta}\_{3,j}\\ does **not** have a \\k\\ subscript.
- \\\boldsymbol{\beta}\_{4,j,k}\mathbf{v}\_i\\: The effect of covariates
  assumed to be **level-specific**. The coefficient vector
  \\\boldsymbol{\beta}\_{4,j,k}\\ **does** have a \\k\\ subscript,
  allowing its effect to be different for each subgroup level. This can
  be seen as an interaction between variable \\x_j\\ and \\v_i\\.

Note that for coefficients \\\boldsymbol{\beta}\_{3,j}\\ and
\\\boldsymbol{\beta}\_{4,j,k}\\ we can use shrinkage or standard priors
depending in our prior beliefs.

#### 3.5.2 Example

Similarly to the example showed for the Global model. The example below
is for the **Age Group** variable, where \\k \in \\ \text{\<65},
\text{≥65} \\\\.

The outcome for a patient *i* in age subgroup *k* is: \\Y\_{ik} \sim
N(\mu\_{ik}, \sigma_k^2)\\.

The linear predictor \\\mu\_{ik}\\ is:

\\ \mu\_{ik} = \underbrace{\beta\_{1,k}}\_{\text{Subgroup Intercept}} +
\underbrace{\beta\_{2,k} \cdot \text{Treatment}\_i}\_{\text{Predictive
Effect}} + \underbrace{\delta_1 \cdot \text{BMI}\_i + \delta_2 \cdot
\text{Smoker}\_i}\_{\text{Constant Prognostic Effects}} +
\underbrace{\delta\_{3,k} \cdot \text{BaselineSBP}\_i}\_{\text{Varying
Prognostic Effect}} \\

- Here, \\\beta\_{1,k}\\ (subgroup-specific intercept) and
  \\\beta\_{2,k}\\ (subgroup-specific predictive effect) are given
  hierarchical priors to borrow information across age groups. For
  example: \\\beta\_{2,k} \sim N(\mu_2, \tau_2^2)\\.
- The effect of Baseline SBP is also allowed to vary by age group
  (\\\delta\_{3,k}\\), while BMI and Smoking effects are assumed
  constant

### 3.6 The Logic of the Intercept

While the intercept (\\\beta_0\\) is not subject to shrinkage (it is a
fixed effect within the `unprogeffect` component), specifying its prior
requires care to ensure the model remains weakly informative without
introducing bias or numerical instability.

**Default Specification:** By default, `bonsaiforest` assigns a weakly
informative Normal prior centered at 0 with a standard deviation of 5:
\\ \beta_0 \sim \mathcal{N}(0, 5^2) \\ This choice is aligned with the
methodology in Wolbers et al. (2025), where a standard deviation of 5 on
the log-hazard or logit scale covers a wide range of plausible effect
sizes (e.g., ratios from 0.08 to 12) while still penalizing extreme,
unrealistic values.

**Linking Prior Scale to Trial Design:** For a more principled choice,
users can adjust the scale parameter (\\\sigma\_{prior}\\) based on
assumptions made during the **trial design and sample size calculation
phase**.

1.  **For Continuous Endpoints:** The trial protocol typically assumes a
    population standard deviation, \\\sigma\_{design}\\, for power
    calculations.
    - **Recommendation:** Set the prior standard deviation to a multiple
      of this design parameter, such as \\\sigma\_{prior} = 10 \times
      \sigma\_{design}\\. This ensures the prior is wide enough to cover
      the plausible range of outcomes envisioned during planning.
2.  **For Count Endpoints (Log Scale):** The intercept represents the
    log of the baseline event rate (assuming the offset is handled
    correctly).
    - **Recommendation:** The default \\\mathcal{N}(0, 5^2)\\ is
      generally appropriate, as it covers a vast multiplicative range of
      event rates. However, if the trial is designed for very rare
      events (e.g., \\\< 0.01\\ events/year) or very frequent events,
      centering the prior on the log of the anticipated event rate from
      the protocol (\\\mu_0 = \log(\lambda\_{design})\\) rather than 0
      can improve convergence.
3.  **For Binary Endpoints (Logit Scale):**
    - **Recommendation:** Use a prior scaled to the **Unit Information
      Standard Deviation (UISD)**, or adhere to the default
      \\\sigma=5\\. A variance of 25 covers the entire range of
      realistic probabilities without concentrating mass on
      deterministic outcomes (0 or 1).

**Special Case: Time-to-Event Endpoints.** For Cox proportional hazards
models, the intercept is **not identifiable** and is **omitted** from
the linear predictor.

- **Mechanism:** The intercept is effectively absorbed into the
  non-parametric baseline hazard function, \\h_0(t)\\.

- **Implementation:** The `brms` backend handles this automatically by
  estimating the baseline hazard (via `bhaz()`) rather than a fixed
  intercept parameter. Therefore, no prior specification for \\\beta_0\\
  is required for survival models.

### 3.7 Endpoint-Specific Models

The package supports four primary endpoint types by specifying the
appropriate likelihood and link function for \\g(\boldsymbol{\mu}) =
\boldsymbol{\eta}\\.

| Endpoint          | Distribution / Model                                       | Link Function                         | Notes                                                                                                    |
|:------------------|:-----------------------------------------------------------|:--------------------------------------|:---------------------------------------------------------------------------------------------------------|
| **Continuous**    | Normal: \\y_i \sim N(\mu_i, \sigma^2)\\                    | Identity: \\\mu_i = \eta_i\\          | \\\sigma\\ can be stratified (e.g., `sigma ~ 0 + clinic_site`).                                          |
| **Binary**        | Bernoulli: \\y_i \sim \text{Bernoulli}(p_i)\\              | Logit: \\\text{logit}(p_i) = \eta_i\\ |                                                                                                          |
| **Count**         | Negative Binomial: \\y_i \sim \text{NegBin}(\mu_i, \phi)\\ | Log: \\\log(\mu_i) = \eta_i\\         | \\\phi\\ (overdispersion) can be stratified. Supports [`offset()`](https://rdrr.io/r/stats/offset.html). |
| **Time-to-event** | Cox Proportional Hazards: \\h_i(t) = h_0(t)\exp(\eta_i)\\  | Log (on the hazard)                   | Intercept is absorbed by the baseline hazard \\h_0(t)\\. \\h_0(t)\\ can be stratified.                   |

### 3.8 Prior Specifications

#### 3.8.1 Default Priors

If priors are not specified, the function applies the following defaults
(as defined in the `fit_brms_model` code):

- **Prognostic Shrunk:** `horseshoe(1)`
- **Prognostic Unshrunk:** `normal(0, 5)`
- **Prognostic Intercept:** `normal(0, 10)`
- **Predictive Shrunk:** `horseshoe(1)`
- **Predictive Unshrunk:** `normal(0, 5)`

**Note:** The intercept prior is *not* applied for
`response_type = "survival"`, as Cox models do not have a global
intercept.

#### 3.8.2 Weakly Informative Priors (Unshrunk Terms)

For all non-penalized coefficients (\\\boldsymbol{\beta\_{1,n}},
\boldsymbol{\beta\_{2,n}}, \boldsymbol{\beta\_{3}}\\), we use weakly
informative priors to aid computational stability without strongly
influencing the posterior.

- **Default:** `normal(0, 10)`. This is a common, weakly regularizing
  prior on the linear predictor scale. Users can specify their own,
  e.g., `student_t(3, 0, 2.5)`.

#### 3.8.3 Regularized Horseshoe Prior (Shrunk Terms)

This is the default shrinkage prior in `bonsaiforest2`, recommended for
its excellent adaptive shrinkage properties ([Wolbers et al.
2025](#ref-wolbers2025using); [Piironen and Vehtari
2017](#ref-piironen2017sparsity)).

- **Concept:** A global-local prior. It has an infinitely tall spike at
  zero (to aggressively shrink noise) and heavy tails (to leave true,
  large signals unshrunk).

- **Hierarchical Specification**:

  \\\begin{aligned} \beta\_{2,l} &\sim N(0, \tau^2 \tilde{\lambda}\_l^2)
  \\ \tilde{\lambda}\_l^2 &= \frac{c^2 \lambda_l^2}{c^2 + \tau^2
  \lambda_l^2} \\ \lambda_l &\sim C^+(0, 1) \quad (\text{Local
  shrinkage}) \\ \tau &\sim C^+(0, \tau_0) \quad (\text{Global
  shrinkage}) \\ c^2 &\sim \text{Inverse-Gamma}(\nu/2, \nu s^2/2)
  \end{aligned}\\

- **Hyperparameter Justification**: The package supports two ways to set
  the crucial global scale \\\tau_0\\:

  1.  **Fixed Default (`scale_global`):** `horseshoe(scale_global = 1)`.
      This is the package default, a robust, general-purpose choice
      ([Wolbers et al. 2025](#ref-wolbers2025using)).
  2.  **Elicited Prior (`par_ratio`):** Sets \\\tau_0\\ based on the
      *prior* belief about the number of effective (non-zero) subgroups,
      \\p\_{eff}\\, out of the total \\L\\ ([Bornkamp
      2025](#ref-bornkamp2025benchmarking); [Piironen and Vehtari
      2017](#ref-piironen2017sparsity)). This is specified via
      \\\text{par\\ratio} = \frac{p\_{eff}}{L - p\_{eff}}\\. This scales
      the prior based on sample size (\\N\\) and the expected sparsity.

**Regularized Horseshoe Implementation in**
[`brms`](https://paulbuerkner.com/brms/)

The [`brms`](https://paulbuerkner.com/brms/) package implements the
Regularized Horseshoe prior using the function
[`horseshoe(...)`](https://paulbuerkner.com/brms/reference/horseshoe.html).
The regularization term (the “slab”) is added to the standard Horseshoe
to ensure robustness and improve computational stability in Stan.

- **Horseshoe Function and Key Arguments**

he Horseshoe prior is set using the function call:

``` r
horseshoe(df = 1, scale_global = 1, df_global = 1, scale_slab = 2, df_slab = 4, par_ratio = NULL, autoscale = TRUE, main = FALSE)
```

| Horseshoe Argument         | Theoretical Equivalent            | Role and Implication                                                                                                                                      |
|:---------------------------|:----------------------------------|:----------------------------------------------------------------------------------------------------------------------------------------------------------|
| \\\mathbf{scale\\global}\\ | \\\tau_0\\                        | Scale of the \\\text{Student-}t\\ prior on the global shrinkage parameter \\\tau\\ Lower values \\\implies\\ **more aggressive global shrinkage**.        |
| \\\mathbf{par\\ratio}\\    | \\\frac{p\_{eff}}{L - p\_{eff}}\\ | **Alternative to** \\\mathbf{scale\\global}\\. Specifies the ratio of expected non-zero to zero coefficients, automatically calculating \\\tau_0\\.       |
| \\\mathbf{df}\\            | DOF for \\\lambda_l\\             | Degrees of freedom of \\\text{Student-}t\\ prior on local shrinkage parameters (\\\lambda_l\\). Default \\\mathbf{df=1}\\ gives the standard half-Cauchy. |
| \\\mathbf{scale\\slab}\\   | \\s\\                             | Scale of the regularization slab. Default \\\mathbf{scale\\slab=2}\\ prevents large coefficients from being over-shrunk.                                  |
| \\\mathbf{autoscale}\\     | Scaling by \\\sigma\\             | Logical. If \\\mathbf{TRUE}\\ (default), the prior is **automatically scaled** using the residual standard deviation (\\\sigma\\).                        |

- **Relationship to Justification:** The \\\mathbf{brms}\\ function
  offers two primary ways to set the global shrinkage level, which
  controls the overall sparsity of the model:

1.  **Fixed Global Scale:** Using the \\\mathbf{scale\\global}\\
    argument, often set to \\\mathbf{1.0}\\ as a general-purpose
    default. This value is internally multiplied by the residual
    standard deviation (\\\sigma\\) when `autoscale = TRUE`.
2.  **Elicited Global Scale (Recommended):** Using the
    \\\mathbf{par\\ratio}\\ argument, which allows the user to specify
    their **prior belief about the sparsity** of the model (i.e., the
    expected number of true non-zero coefficients, \\p\_{eff}\\). This
    is the recommended approach as it accounts for the number of
    coefficients \\L\\ and the sample size \\N\\ to determine the
    optimal \\\tau_0\\.

- **Choosing** \\\mathbf{scale\\global}\\ **and**
  \\\mathbf{par\\ratio}\\: The choice of \\\tau_0\\ (controlled by
  \\\mathbf{scale\\global}\\ or \\\mathbf{par\\ratio}\\) is crucial, as
  it sets the overall magnitude of shrinkage.

  - **Fixed Default (**\\\mathbf{scale\\global=1}\\):

    - **Use when:** You want a simple, robust, general-purpose prior
      without explicit sparsity knowledge.
    - **Caveat:** May result in insufficient shrinkage if the true
      number of non-zero coefficients is small.

  - **Elicited Sparsity (**\\\mathbf{par\\ratio}\\):

    - **Calculation:** \\\text{par\\ratio} = \frac{p\_{eff}}{L -
      p\_{eff}}\\, where \\p\_{eff}\\ is the expected number of non-zero
      coefficients and \\L\\ is the total number of coefficients.
    - **Use when:** You have a strong prior belief about the number of
      relevant effects. This is generally preferred for adaptive
      shrinkage as it scales \\\tau_0\\ appropriately.

|  \\\mathbf{par\\ratio}\\ Value  | Expected \\p\_{eff}\\ (Sparsity) | Interpretation of Shrinkage                                                                          |
|:-------------------------------:|:---------------------------------|:-----------------------------------------------------------------------------------------------------|
|     **Small (e.g., 0.05)**      | Very Few Non-Zero Coefficients   | **High Shrinkage:** \\\tau_0\\ is small, aggressively shrinking coefficients towards zero.           |
|    **Moderate (e.g., 0.2)**     | A Moderate Fraction is Non-Zero  | **Moderate Shrinkage:** A less aggressive \\\tau_0\\, allowing more coefficients to remain unshrunk. |
| \\\mathbf{par\\ratio}\\ Ignored | \\\mathbf{scale\\global}\\ used  | Uses a fixed \\\tau_0\\ that doesn’t explicitly depend on expected sparsity.                         |

> **Recommendation:** While \\\mathbf{scale\\global=1}\\ is the quick
> default, specifying \\\mathbf{par\\ratio}\\ based on expert knowledge
> of expected sparsity is the **most theoretically sound and adaptive**
> method for setting the global shrinkage scale.

#### 3.8.4 R2D2 Prior (Shrunk Terms)

- **Concept:** This prior is uniquely derived by first placing a prior
  on the model’s coefficient of determination, \\R^2\\, which quantifies
  the proportion of variance explained by the predictors. This is
  arguably more intuitive than specifying priors directly on
  coefficients.

- **Hierarchical Specification (Marginal Version):** \\ \begin{aligned}
  \beta\_{2,l} &\sim DE(\sigma \sqrt{\phi_l \omega / 2}) \quad
  (\text{Double-Exponential/Laplace kernel}) \\ \boldsymbol{\phi} &=
  (\phi_1, \dots, \phi_L) \sim \text{Dirichlet}(a\_\pi, \dots, a\_\pi)
  \quad (\text{Local shrinkage}) \\ \omega &\sim \text{Beta-Prime}(a, b)
  \quad (\text{Global shrinkage}) \end{aligned} \\ This is equivalent to
  assuming that the proportion of variance explained by the interaction
  terms follows \\R^2 = \frac{\omega}{1+\omega} \sim \text{Beta}(a,b)\\.

- **Justification of Hyperparameters:**

  - Zhang et al. ([2022](#ref-zhang2022bayesian)) provide clear
    guidance. A fully automatic approach is to set \\b=0.5\\ to achieve
    **Cauchy-like heavy tails**.

  - The concentration parameter \\a\_\pi\\ controls sparsity. A smaller
    \\a\_\pi\\ (e.g., 0.2) implies **higher shrinkage**, concentrating
    the prior variance on fewer coefficients, while a larger \\a\_\pi\\
    (e.g., 0.5) spreads the variance more evenly, implying **lower
    shrinkage**.

  - **Default Recommendation:** Set \\b=0.5\\ and offer options for the
    user based on desired shrinkage strength:

    - **High Shrinkage:** \\a\_\pi=0.2\\
    - **Low Shrinkage:** \\a\_\pi=0.5\\

**R2D2 Prior Implementation in**
[`brms`](https://paulbuerkner.com/brms/)

The [`brms`](https://paulbuerkner.com/brms/) package implements this
prior using the function
[`R2D2(...)`](https://paulbuerkner.com/brms/reference/R2D2.html),
providing a high-level parametrization that is more intuitive than
setting the parameters \\a\\, \\b\\, and \\a\_\pi\\ directly. This
function translates desired \\\mathbf{R^2}\\ properties and shrinkage
strength into the underlying prior distribution’s parameters.

- **R2D2 Function and Key Arguments**

The R2D2 prior is set using the function call:

``` r
R2D2(mean_R2 = 0.5, prec_R2 = 2, cons_D2 = 0.5, autoscale = TRUE, main = FALSE)
```

| R2D2 Argument          | Theoretical Equivalent  | Role and Implication                                                                                                                                                          |
|:-----------------------|:------------------------|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| \\\mathbf{mean\\R2}\\  | \\\mathbb{E}\[R^2\]\\   | Sets the **prior expected value** for the model’s fixed effect \\R^2\\.                                                                                                       |
| \\\mathbf{prec\\R2}\\  | \\a+b\\ (Concentration) | Controls the **precision** of the Beta prior on \\R^2\\. Lower values mean a flatter prior.                                                                                   |
| \\\mathbf{cons\\D2}\\  | \\a\_{\pi}\\            | The **concentration parameter** for the \\\text{Dirichlet}\\ distribution on the variance components (local shrinkage). **Lower values imply higher shrinkage/sparsity.**     |
| \\\mathbf{autoscale}\\ | Scaling by \\\sigma\\   | Logical. If \\\mathbf{TRUE}\\ (default), the prior is **automatically scaled** using the residual standard deviation (\\\sigma\\). This ensures the prior is scale-invariant. |

- **Relationship to Justification:** The
  [`brms`](https://paulbuerkner.com/brms/) function parameters are
  directly related to the justified hierarchical parameters:

  - The \\\mathbf{cons\\D2}\\ parameter is the direct counterpart to the
    sparsity parameter \\\mathbf{a\_{\pi}}\\. The default recommendation
    \\\mathbf{cons\\D2 = 0.5}\\ corresponds precisely to the
    \\\mathbf{a\_{\pi}=0.5}\\ low shrinkagesetting.
  - The \\\mathbf{mean\\R2}\\ and \\\mathbf{prec\\R2}\\ arguments define
    the \\\text{Beta}(a,b)\\ prior on \\R^2\\, where: \\a =
    \text{mean_R2} \times \text{prec_R2}\\ \\b = (1 - \text{mean_R2})
    \times \text{prec_R2}\\ The choice of \\\mathbf{prec\\R2=2}\\ (the
    default) alongside \\\mathbf{mean\\R2=0.5}\\ results in \\a=1\\ and
    \\b=1\\, which is a **Uniform prior** on \\R^2 \sim \text{Beta}(1,
    1)\\, representing a maximally weak prior belief on global fit.
  - The essential \\\mathbf{b=0.5}\\ setting for **Cauchy-like heavy
    tails** (preventing over-shrinkage) is **implicitly** handled by the
    R2D2 function’s underlying structure and does not require manual
    specification.

- **Choosing** \\\mathbf{mean\\R2}\\ and \\\mathbf{prec\\R2}\\ :
  Choosing these parameters involves quantifying your prior knowledge
  about the model’s overall predictive power before observing the data.
  They define a \\\text{Beta}(a,b)\\ distribution on \\R^2\\.

  - Choosing \\\mathbf{mean\\R2}\\ (Prior Expected \\R^2\\): This sets
    the **center** of your prior belief.

    - **Default (0.5):** Use if you have no strong prior knowledge. It
      implies a neutral expectation that the fixed effects will explain
      half of the variance.

    - **Informative Low (e.g., 0.1 to 0.3):** Use if you expect a very
      weak or noisy relationship.

    - **Informative High (e.g., 0.7 to 0.9):** Use if you expect a
      highly predictive model based on strong prior evidence.

  - Choosing \\\mathbf{prec\\R2}\\ (Prior Confidence in \\R^2\\): This
    determines the **concentration** or **certainty** of the Beta prior
    around the chosen \\\text{mean_R2}\\.

[TABLE]

> **Recommendation:** The **Default** \\\mathbf{mean\\R2 = 0.5}\\ and
> \\\mathbf{prec\\R2 = 2}\\ provides the most flexible, weakly
> informative setting for \\R^2\\, which is generally recommended unless
> strong prior knowledge dictates otherwise.

#### 3.8.5 Normal Hierarchical Prior (for One-Variable-at-a-Time Model)

This is the standard hierarchical model for subgroup analysis, which
assumes that treatment effects for levels within a subgrouping variable
are drawn from a common normal distribution. It’s less complex than the
Horseshoe or R2D2 priors but provides effective shrinkage, especially
when the number of subgroup levels is small.

- **Hierarchical Specification:** For a given subgrouping variable \\j\\
  with levels \\k=1, \dots, K_j\\:

\\\begin{aligned} \beta\_{2,j,k} &\sim \mathcal{N}(\mu\_{2,j},
\tau^2\_{2,j}) \quad (\text{Level-specific treatment effects}) \\
\mu\_{2,j} &\sim \mathcal{N}(0, s^2) \quad (\text{Common mean effect})
\\ \tau\_{2,j} &\sim \text{Half-Normal}(a) \quad (\text{Between-subgroup
heterogeneity}) \end{aligned} \\

- **Justification of Hyperparameters:** The key is setting the scale
  parameter a for the Half-Normal prior on the between-subgroup standard
  deviation, \\\tau\_{2,j}\\, as this controls the degree of shrinkage.
  As recommended by Bornkamp ([2025](#ref-bornkamp2025benchmarking)),
  this can be linked to the planned treatment effect (
  \\\delta\_{plan}\\) from the trial protocol, making the prior choice
  interpretable and consistent across different endpoints. The choice of
  a implies a prior on the plausible difference in treatment effects
  between any two subgroups.

### 3.9 Estimation (MCMC)

The joint posterior distribution of all parameters is complex and has no
closed-form solution. We use Markov Chain Monte Carlo (MCMC) methods to
draw samples from this distribution. Specifically, the package uses the
**No-U-Turn Sampler (NUTS)**, an efficient variant of Hamiltonian Monte
Carlo (HMC), as implemented in **Stan** via the `brms` package. The
output is a set of posterior samples (e.g., 4000 draws) representing the
joint posterior distribution.

### 3.10 Standardization for Marginal Effects (G-computation)

To obtain interpretable *marginal* effects for each subgroup, the
package implements a standardization (G-computation) procedure. This
process is repeated for *each* MCMC draw to generate a full posterior
distribution of the marginal effect.

**For a single MCMC draw, the step-by-step process is:**

1.  **Select Parameters:** Take one draw from the joint posterior
    distribution for all model parameters (all \\\boldsymbol{\beta}\\’s,
    \\\sigma\\, \\\phi\\, etc.).

2.  **Create Counterfactuals:** For *every patient* \\i\\ in the
    original dataset, create two scenarios:

    - Scenario A: Patient \\i\\’s covariates, with treatment \\s_i=0\\
      (Control).
    - Scenario B: Patient \\i\\’s covariates, with treatment \\s_i=1\\
      (Treatment).

    For Scenario A (control), all treatment-by-subgroup interaction
    dummies are set to 0. For Scenario B (treatment), each interaction
    dummy `trt_subgroupLEVEL` is set to 1 if the patient belongs to that
    subgroup level, 0 otherwise.

3.  **Predict Outcomes:** Use the model formula and the parameters from
    Step 1 to predict the outcome for *every patient* under both
    Scenario A (\\\hat{\mu}\_{i,0}\\) and Scenario B
    (\\\hat{\mu}\_{i,1}\\). The
    [`brms::posterior_epred()`](https://mc-stan.org/rstantools/reference/posterior_epred.html)
    function automatically handles the computation using the model’s
    design matrix and parameters.

4.  **Average within Subgroups:** For a subgroup of interest \\l\\
    (e.g., “Female”):

    - Calculate the average predicted outcome under control:
      \\\hat{\mu}\_{l,0} = \text{mean}(\hat{\mu}\_{i,0})\\ for all \\i\\
      where \\x\_{il}=1\\.
    - Calculate the average predicted outcome under treatment:
      \\\hat{\mu}\_{l,1} = \text{mean}(\hat{\mu}\_{i,1})\\ for all \\i\\
      where \\x\_{il}=1\\.

5.  **Calculate Effect Measure:** Compute the marginal effect for
    subgroup \\l\\ from these averaged predictions. This depends on the
    endpoint type:

    - **Continuous:** Mean Difference = \\\hat{\mu}\_{l,1} -
      \hat{\mu}\_{l,0}\\.
    - **Binary:** Odds Ratio = \\\frac{\hat{\mu}\_{l,1} / (1 -
      \hat{\mu}\_{l,1})}{\hat{\mu}\_{l,0} / (1 - \hat{\mu}\_{l,0})}\\
      (where \\\hat{\mu}\\ is the predicted probability).
    - **Count:** Rate Ratio = \\\hat{\mu}\_{l,1} / \hat{\mu}\_{l,0}\\
      (where \\\hat{\mu}\\ is the predicted rate).
    - **Time-to-Event:** Average Hazard Ratio (AHR). This is the most
      complex, as it requires predicting the full survival curve
      \\S_i(t)\\ for each patient. The marginal survival curves
      \\\hat{S}\_{l,0}(t)\\ and \\\hat{S}\_{l,1}(t)\\ are calculated,
      and the AHR is computed from them, as the marginal curves may not
      be proportional.

6.  **Repeat:** This process (Steps 1-5) is repeated for all MCMC draws,
    yielding a full posterior distribution (e.g., 4000 values) for the
    marginal effect of each subgroup. The final point estimate (median)
    and 95% credible interval are taken from this distribution.

## 4 Mapping of Statistical Methods to `bonsaiforest2` Functions

This section connects the statistical methodology (Section 3) to the
package functions. The package provides a three-level architecture:
high-level user interface functions, modular worker functions for
advanced control, and internal helpers.

### 4.1 High-Level User Interface

These are the main functions users typically call in their analysis
workflow:

- **[`run_brms_analysis()`](https://openpharma.github.io/bonsaiforest2/reference/run_brms_analysis.md)**

  - **Maps to:** Sections 3.2 (Models), 3.3 (Endpoints), 3.4 (Priors),
    3.5 (Estimation).
  - **Action:** Main entry point that orchestrates the complete modeling
    workflow.
    - Internally calls
      [`prepare_formula_model()`](https://openpharma.github.io/bonsaiforest2/reference/prepare_formula_model.md)
      to construct the three-component `brmsformula` structure and
      validate inputs.
    - Then calls
      [`fit_brms_model()`](https://openpharma.github.io/bonsaiforest2/reference/fit_brms_model.md)
      to run MCMC sampling via
      [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html)
      and Stan.
    - Stores essential metadata (treatment variable, response type,
      original data) as attributes on the fitted model object for
      downstream analysis.
    - Handles all formula components: `unshrunktermeffect` (unshrunk
      terms), `shprogeffect` (shrunk prognostic), and `shpredeffect`
      (shrunk predictive).
    - Returns a `brmsfit` object with enriched attributes.

- **[`summary_subgroup_effects()`](https://openpharma.github.io/bonsaiforest2/reference/summary_subgroup_effects.md)**

  - **Maps to:** Section 3.6 (Standardization / G-computation).
  - **Action:** User-facing wrapper for marginal effect summarization.
    - Internally calls
      [`estimate_subgroup_effects()`](https://openpharma.github.io/bonsaiforest2/reference/estimate_subgroup_effects.md)
      to perform the G-computation procedure.
    - Automatically extracts treatment variable, response type, and data
      from model attributes (set by
      [`fit_brms_model()`](https://openpharma.github.io/bonsaiforest2/reference/fit_brms_model.md)).
    - Auto-detects subgroups from all formula components by searching
      for treatment interactions in `unshrunktermeffect`,
      `shprogeffect`, and `shpredeffect`.
    - Returns a `subgroup_summary` S3 object containing estimates,
      credible intervals, and metadata.

- **[`plot()`](https://rdrr.io/r/graphics/plot.default.html)**

  - **Maps to:** Final visualization of Section 3.6 results.
  - **Action:** S3 method for `subgroup_summary` objects. Generates
    publication-ready forest plots from marginal effect estimates with
    median and 95% credible intervals.

### 4.2 Advanced/Modular Interface

For users requiring fine-grained control over the modeling process, the
package exposes intermediate worker functions:

- **[`prepare_formula_model()`](https://openpharma.github.io/bonsaiforest2/reference/prepare_formula_model.md)**

  - **Parameters:**
    - `data`: Dataset containing all variables
    - `response_formula`: Response specification (e.g., `outcome ~ trt`
      or `Surv(time, status) ~ trt`)
    - `unshrunk_terms_formula`: Unshrunk terms specification (may
      include main effects and treatment interactions)
    - `shrunk_prognostic_formula`: Prognostic main effects to be
      regularized (optional)
    - `shrunk_predictive_formula`: Treatment interactions to be
      regularized (optional)
    - `response_type`: One of `"binary"`, `"count"`, `"continuous"`, or
      `"survival"`
    - `stratification_formula`: Stratification variable (optional)
  - **Action:** Constructs and validates the three-component
    `brmsformula` structure:
    - `unshrunktermeffect`: Unshrunk terms with weakly informative
      priors (intercept, main treatment effect, pre-specified
      covariates). Specified via `unshrunk_terms_formula`.
    - `shprogeffect`: Shrunk prognostic effects with strong
      regularization priors. Specified via `shrunk_prognostic_formula`.
    - `shpredeffect`: Shrunk predictive effects (treatment interactions)
      with strong regularization priors. Specified via
      `shrunk_predictive_formula`.
  - Supports flexible interaction syntax: colon notation
    (`trt:subgroup`), star notation (`trt*subgroup`), or random effects
    notation (`(trt || subgroup)`), which can be mixed in the same
    model.
  - Validates model hierarchy: checks if predictive terms have
    corresponding prognostic main effects, issuing warnings if
    violations are detected (star notation automatically includes main
    effects and is excluded from this check).
  - Applies automatic contrast coding: one-hot encoding for shrunk terms
    (all levels represented for exchangeability) and dummy encoding for
    unshrunk terms (reference level dropped).
  - Constructs design matrices \\\mathbf{X_n}, \mathbf{X_p},
    \mathbf{Z_n}, \mathbf{Z_p}, \mathbf{U}\\ from formula
    specifications.
  - Returns a list containing the complete `brmsformula`, transformed
    data with proper contrasts, response type, treatment variable name,
    Stan parameter names, and metadata.

- **[`fit_brms_model()`](https://openpharma.github.io/bonsaiforest2/reference/fit_brms_model.md)**

  - **Action:** Executes the MCMC sampling and model fitting.
    - Takes output from
      [`prepare_formula_model()`](https://openpharma.github.io/bonsaiforest2/reference/prepare_formula_model.md)
      or manually specified formula components.
    - Assigns priors to each formula component based on user
      specifications or defaults.
    - Calls
      [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html)
      to compile Stan code and run MCMC chains.
    - Handles endpoint-specific configurations (Section 3.3): likelihood
      functions, link functions, and stratification.
    - Attaches essential metadata as attributes to the fitted model for
      use by
      [`estimate_subgroup_effects()`](https://openpharma.github.io/bonsaiforest2/reference/estimate_subgroup_effects.md).

- **[`estimate_subgroup_effects()`](https://openpharma.github.io/bonsaiforest2/reference/estimate_subgroup_effects.md)**

  - **Action:** Implements the complete G-computation procedure
    described in Section 3.8:
    1.  For each MCMC draw, creates counterfactual datasets (all
        patients as treated vs. control).
    2.  Generates posterior predictions using
        [`brms::posterior_epred()`](https://mc-stan.org/rstantools/reference/posterior_epred.html)
        or `brms::posterior_survfit()`.
    3.  Averages predictions within each subgroup to obtain marginal
        expected outcomes.
    4.  Computes endpoint-specific effect measures:
        - Continuous: Mean difference (\\\hat{\mu}\_{l,1} -
          \hat{\mu}\_{l,0}\\)
        - Binary: Odds ratio
        - Count: Rate ratio
        - Survival: Average hazard ratio (AHR)
  - Auto-detects subgroups by parsing all formula components for
    treatment interactions (both colon `:` syntax and pipe `||` syntax
    for random effects).
  - Returns raw posterior draws and summary statistics (median, credible
    intervals).

### 4.3 Internal Helper Functions

These functions support the main workflow but are not typically called
directly by users:

- **`.prepare_subgroup_vars()`**: Detects subgroup variables from model
  formula by identifying treatment interaction terms.
- **`.create_counterfactual_data()`**: Generates counterfactual datasets
  for G-computation.
- **`.compute_effect_measure()`**: Calculates appropriate effect measure
  based on response type.

## References

Bornkamp, Björn. 2025. “Benchmarking Bayesian subgroup shrinkage methods
on clinical data.” In *PSI Conference 2025*. London, UK.

Piironen, J., and A. Vehtari. 2017. “Sparsity information and
regularization in the horseshoe and other shrinkage priors.” *Electronic
Journal of Statistics* 11: 5018–51.

Wang, Yun, Wenda Tu, William Koh, James Travis, Robert Abugov, Kiya
Hamilton, Mengjie Zheng, Roberto Crackel, Pablo Bonangelino, and Mark
Rothmann. 2024. “Bayesian hierarchical models for subgroup analysis.”
*Pharmaceutical Statistics* 23: 1065–83.

Wolbers, Marcel, Mar Vázquez Rabuñal, Ke Li, Kaspar Rufibach, and Daniel
Sabanés Bové. 2025. “Using shrinkage methods to estimate treatment
effects in overlapping subgroups in randomized clinical trials with a
time-to-event endpoint.” *Statistical Methods in Medical Research*,
1–12.

Zhang, Yan Dora, Brian P. Naughton, Howard D. Bondell, and Brian J.
Reich. 2022. “Bayesian regression using a prior on the model fit: The
R2-D2 shrinkage prior.” *Journal of the American Statistical
Association* 117 (538): 862–74.
