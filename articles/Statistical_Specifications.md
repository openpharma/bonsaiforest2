# Statistical Specifications

## 1. Scope of this document

This document describes the statistical methods implemented in the
`bonsaiforest2` R package for Bayesian subgroup analysis of continuous,
binary, count, and time-to-event outcomes in randomized clinical trials
(RCTs).

The package unabled the implementation of both a **global modeling
approach** (Wolbers et al. 2025), that estimates all main (prognostic)
and interaction (predictive) effects in a single model, and a
**One-variable-at a time** (Wang et al. 2024), that estimates the
interaction effects using a different model for each subgrouping
variable. This contrasts with traditional methods that analyze subgroups
separately or one variable at a time.

The core features described in this document are:

1.  **Flexible Model Specification:** A unified framework to distinguish
    between:
    - **Prognostic** (main) vs. **Predictive** (interaction) effects.
    - **Shrunk** vs. **Unshrunk** coefficients, allowing exploratory
      terms to be penalized while pre-specified terms are not.
2.  **Advanced Shrinkage Priors:** Implementation of state-of-the-art
    shrinkage priors, including the **Regularized Horseshoe** (Piironen
    and Vehtari 2017; Wolbers et al. 2025) and **R2D2** (Zhang et al.
    2022), to provide robust control of false-positive findings while
    identifying potential heterogeneity.
3.  **Marginal Effect Estimation:** Use of **standardization
    (G-computation)** to derive interpretable *marginal* subgroup
    treatment effects, which average over the covariate distributions
    within each subgroup (Wolbers et al. 2025).

This document is structured as follows: Section 2 provides a conceptual
introduction to subgroup analysis challenges and the Bayesian shrinkage
solution. Section 3 details the statistical methodology, including
notation, the global and the one-variable at a time models,
endpoint-specific likelihoods, prior specifications, and the
standardization procedure. Section 4 maps these statistical methods to
the core functions in the `bonsaiforest2` package.

## 2. Introduction to Subgroup Analysis and Estimation

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

`bonsaiforest2` solves this by using **Standardization
(G-computation)**. This procedure correctly averages over the specific
covariate distribution of each subgroup to estimate a valid marginal
treatment effect, ensuring all subgroups are compared on the same
interpretable scale.

## 3. Statistical Methodology

### 3.1 Setting and Notation

We establish the following notation for the global model:

- $i = 1,\ldots,N$: Index for individual patients.
- $y_{i}$: The outcome for patient $i$.
- $s_{i}$: The binary treatment indicator ($s_{i} = 1$ for treatment,
  $s_{i} = 0$ for control).
- $l = 1,\ldots,L$: Index for an overall subgroup level (e.g., “Male”,
  “Female”, “EU”, “NA”).
- $x_{il}$: Indicator (0/1) for patient $i$ belonging to subgroup level
  $l$ (used for prognostic effects).
- $z_{il} = s_{i} \cdot x_{il}$: The interaction between treatment and
  subgroup level $l$ (used for predictive effects).
- $u_{iv}$: Value of the $v$-th additional non-subgroup covariate for
  patient $i$.
- $\mathbf{μ}$: The $N \times 1$ vector of expected outcomes.
- $g( \cdot )$: A link function (e.g., logit, log).
- $\mathbf{η}$: The $N \times 1$ vector of the linear predictor.

### 3.2. Design Matrix Considerations

A critical implementation detail is the construction of the design
matrix for the subgroup-by-treatment interactions.

We will then follow the same idea as in the `bonsaiforest` library. We
will treat prognostic effects and predictive effects differently, as
they serve different purposes in the model:

- **For prognostic effects :** The goal is interpretability and avoiding
  collinearity with the intercept. Here, standard dummy coding is the
  appropriate choice. For each of the J subgrouping variables, one level
  is dropped as a reference, resulting in L−J total columns for the
  subgroup main effects.

- **For predictive effects:** The goal is to treat all subgroup levels
  symmetrically under the assumption of exchangeability. Using one-hot
  encoding is essential here. This ensures that the shrinkage prior
  pulls all predictive effects toward a common center without
  privileging a reference level.

### 3.3 The Global Model

The package enables the implementation of a global model (Wolbers et al.
2025) where all effects are estimated in a single, comprehensive
regression formula.

#### 3.3.1 Full Model Equation

***Matrix form:***

The model for the linear predictor, $\mathbf{η}$, for all $N$ patients
is given by:

$$g({\mathbf{μ}}) = {\mathbf{η}} = \underset{\text{Penalized prognostic effects}}{\underbrace{\mathbf{X}_{\mathbf{p}}{\mathbf{β}}_{\mathbf{1},\mathbf{p}}}} + \underset{\text{Non-penalized prognostic effects}}{\underbrace{\mathbf{X}_{\mathbf{n}}{\mathbf{β}}_{\mathbf{1},\mathbf{n}}}} + \underset{\text{Penalized predictive effects}}{\underbrace{\mathbf{Z}_{\mathbf{p}}{\mathbf{β}}_{\mathbf{2},\mathbf{p}}}} +$$$$\underset{\text{Non-penalized predictive effects}}{\underbrace{\mathbf{Z}_{\mathbf{n}}{\mathbf{β}}_{\mathbf{2},\mathbf{n}}}} + \underset{\text{Extra prognostic effects}}{\underbrace{\mathbf{U}{\mathbf{β}}_{\mathbf{3}}}}$$

The Bayesian hierarchical structure is applied as follows:

- **Penalized Coefficients**: ${\mathbf{β}}_{\mathbf{1},\mathbf{p}}$ and
  ${\mathbf{β}}_{\mathbf{2},\mathbf{p}}$ are given shrinkage priors.

- **Non-Penalized Coefficients**: ${\mathbf{β}}_{\mathbf{1},\mathbf{n}}$
  and ${\mathbf{β}}_{\mathbf{2},\mathbf{n}}$ are given standard, weakly
  informative priors.

- **Extra prognostic effects Coefficients**: ${\mathbf{β}}_{\mathbf{3}}$
  are given standard, weakly informative priors.

**Components definition:**

- $\mathbf{μ}$ is the $N \times 1$ vector of expected outcomes for all
  patients.

- $g( \cdot )$ is the link function (e.g., logit for binary data).

- $\mathbf{X} = \left\lbrack \mathbf{X}_{\mathbf{p}}\;\;\mathbf{X}_{n} \right\rbrack$
  is the $N \times M$ design matrix for the **prognostic
  effects**[¹](#fn1). It typically includes:

  - An intercept column. (I guess we won’t like to penalize the
    intercept so it will be contained in $\mathbf{X}_{n}$).

  - A column for the main treatment effect.

  - Columns for the prognostic effects of each subgroup level (using
    dummy coding, so one reference level per factor is dropped). Each
    prognostic effect will be contained in $\mathbf{X}_{n}$ or
    $\mathbf{X}_{p}$ depending if we want to shrink or not the effect.

    Example:

  | Patient | Intercept | trt | sexF | regionUS | regionAsia |
  |:--------|:---------:|:---:|:----:|:--------:|:-----------|
  | 1       |     1     |  0  |  0   |    0     | 0          |
  | 2       |     1     |  1  |  1   |    0     | 1          |
  | 3       |     1     |  0  |  1   |    1     | 0          |
  | 4       |     1     |  1  |  0   |    1     | 0          |

  An example design matrix for a global subgroup model.

  *Note:* M=1+1+L-J= 1 degree for the intercept+ 1 for the treatment
  effect + L (Number of subgroups) - J(Number of variables) because we
  get the reference levels out.

- ${\mathbf{β}}_{\mathbf{1},\mathbf{p}}$ is the $P \times 1$ vector of
  coefficients for the prognostic effects that we want to shrink. We
  will give an shrinkage prior to this coefficients.

- ${\mathbf{β}}_{\mathbf{1},\mathbf{n}}$ is the $(M - P) \times 1$
  vector of coefficients for the prognostic effects that we don’t want
  to shrink.

- $\mathbf{Z} = \left\lbrack \mathbf{Z}_{\mathbf{p}}\;\;\mathbf{Z}_{n} \right\rbrack$
  is the $N \times L$ design matrix for the **predictive effects**. Its
  columns represent the interaction between the treatment and each of
  the $L$ subgroup levels across all subgrouping variables. This matrix
  is constructed without dropping a reference level to treat all
  subgroups symmetrically.

| Patient | trt x sexM | trt x sexF | trt x regionEU | trt x regionUS | trt x regionAsia |
|:-------:|:----------:|:----------:|:--------------:|:--------------:|:----------------:|
|    1    |     0      |     0      |       0        |       0        |        0         |
|    2    |     0      |     1      |       1        |       0        |        0         |
|    3    |     0      |     0      |       0        |       0        |        0         |
|    4    |     1      |     0      |       0        |       1        |        0         |
|    5    |     0      |     0      |       0        |       0        |        0         |
|    6    |     0      |     1      |       0        |       0        |        1         |

- ${\mathbf{β}}_{\mathbf{2},\mathbf{p}}$ is the $R \times 1$ vector of
  coefficients for the predictive effects that we want to shrink. We
  will give an shrinkage prior to this coefficients.
- ${\mathbf{β}}_{\mathbf{2},\mathbf{n}}$ is the $(L - R) \times 1$
  vector of coefficients for the predictive effects that we don’t want
  to shrink.
- $U$ is the $N \times V$ design matrix for the extra variables we want
  to take into account in the fitting of the model. This will also be
  **unpenalized.**
- ${\mathbf{β}}_{\mathbf{3}}$ is the $V \times 1$ vector of extra
  prognostic effects coefficients

|                       | **Shrinkage**                                                 | **Not Shrinkage**                                             |
|:----------------------|:--------------------------------------------------------------|:--------------------------------------------------------------|
| **Prognostic Effect** | $\mathbf{X}_{\mathbf{p}}{\mathbf{β}}_{\mathbf{1},\mathbf{p}}$ | $\mathbf{X}_{\mathbf{n}}{\mathbf{β}}_{\mathbf{1},\mathbf{n}}$ |
| **Predictive Effect** | $\mathbf{Z}_{\mathbf{p}}{\mathbf{β}}_{\mathbf{2},\mathbf{p}}$ | $\mathbf{Z}_{\mathbf{n}}{\mathbf{β}}_{\mathbf{2},\mathbf{n}}$ |

***Linear form:***

To make it more concrete, the linear predictor for a single patient $i$
can be written as:

$$\eta_{i} = \underset{\text{Prognostic Effects}}{\underbrace{\beta_{0} + \beta_{\text{treat}}s_{i} + \underset{\text{penalized}}{\underbrace{\sum\limits_{l = 1}^{P}\beta_{1,p}^{l}x_{n}^{il}}} + \underset{\text{Not penalized}}{\underbrace{\sum\limits_{l = P + 1}^{M}\beta_{1,n}^{l}x_{n}^{il}}} + \underset{\text{Extra}\mspace{6mu}}{\underbrace{\sum\limits_{v = 1}^{V}\beta_{3}^{v}u^{iv}}}}} + \underset{\text{Predictive effects}}{\underbrace{\underset{\text{penalized}\mspace{6mu}}{\underbrace{\sum\limits_{l = 1}^{L - R}\beta_{2,p}^{l}z_{p}^{il}}} + \underset{\text{Not penalized}}{\underbrace{\sum\limits_{l = L - R + 1}^{L}\beta_{2,n}^{l}z_{n}^{il}}}}}$$

**Note**: See that even if $z^{il} = x^{il} \cdot s^{i}$, as the
prognostic effect or the predicitve effect of the same subgrouping
variable may be penalized or not, we separate it in $X$ and $z$ in the
formula.

#### 3.3.2 Example

A randomized, double-blind trial comparing a new drug, **DrugX**,
against a **placebo** in patients with hypertension.

- **Primary Endpoint (**$Y_{i}$): Change from baseline in systolic blood
  pressure (SBP) in mmHg at 6 months.
- **Subgroups of Interest**:
  - **Age Group**: `<65` vs. `≥65` years.
  - **Sex**: `Male` vs. `Female`.
  - **Kidney Function**: eGFR `<60` (impaired) vs. `≥60` (normal).
- **Covariates**: `BaselineSBPᵢ`, `BMIᵢ`, `Smokerᵢ` (binary).

**Global Model with Shrinkage**

The outcome is assumed to be normally distributed:
$Y_{i} \sim N\left( \mu_{i},\sigma^{2} \right)$.

The linear predictor $\mu_{i}$ is modeled as:

$$\mu_{i} = \underset{\text{Intercept}}{\underbrace{\beta_{0}}} + \underset{\text{Main Trt Effect}}{\underbrace{\beta_{\text{treat}} \cdot \text{Treatment}_{i}}} + \underset{\text{Prognostic Subgroup Effects}}{\underbrace{\sum\limits_{k = 1}^{K}\alpha_{k} \cdot \text{Subgroup}_{ik}}} + \underset{\text{Shrunken Predictive Effects}}{\underbrace{\sum\limits_{k = 1}^{K}\gamma_{k} \cdot \left( \text{Subgroup}_{ik} \cdot \text{Treatment}_{i} \right)}} + \underset{\text{Extra Prognostic effects}}{\underbrace{\sum\limits_{j = 1}^{J}\delta_{j} \cdot \text{Covariate}_{ij}}}$$

$$= \underset{\text{Intercept}}{\underbrace{\beta_{0}}} + \underset{\text{Main Trt Effect}}{\underbrace{\beta_{\text{treat}} \cdot \text{Treatment}_{i}}} + \underset{\text{Prognostic Subgroup Effects}}{\underbrace{\alpha_{1} \cdot \text{Age Group}_{i} + \alpha_{2} \cdot \text{Sex}_{i} + \alpha_{3} \cdot \text{eGFR}_{i}}} +$$

$$+ \underset{\text{Predictive Effects}}{\underbrace{\gamma_{1} \cdot \left( \text{Age Group}_{i} \cdot \text{Treatment}_{i} \right) + \gamma_{2} \cdot \left( \text{Sex}_{i} \cdot \text{Treatment}_{i} \right) + \gamma_{3} \cdot \left( \text{eGFR}_{i} \cdot \text{Treatment}_{i} \right)}} + \underset{\text{Extra Prognostic effects}}{\underbrace{\delta_{1} \cdot \text{BaselineSBP}_{i} + \delta_{2} \cdot \text{BMI}_{i} + \delta_{3} \cdot \text{Smoker}_{i}}}$$

### 3.4. One-Variable-at-a-Time Model

The **one-variable-at-a-time model** (Wang et al. 2024) is a Bayesian
hierarchical model (BHM) that improves upon standard subgroup analysis
by borrowing information across subgroup levels to produce more stable
and reliable estimates. Instead of analyzing each subgroup level in
isolation, it analyzes all levels *within a single subgrouping variable*
(e.g., all age groups) together in one model.

For example, for the variable *Sex* (with levels Male and Female), you
would fit one hierarchical model. In this single model, the treatment
effects for Male ($\beta_{2,\text{sex,male}}$) and Female
($\beta_{2,\text{sex,female}}$) are treated as exchangeable and are
linked by a common prior distribution. Then, a completely separate
hierarchical model would be fit for the variable *Region*.

#### 3.4.1 Full Model Equation

For a given subgrouping variable $j$, we can specify a flexible linear
predictor that allows for some covariate effects to be constant across
the levels of $j$, while others are allowed to vary.

Let’s partition the patient-level covariates into two sets:

- $\mathbf{u}_{i}$: A vector of covariates assumed to have a **constant
  prognostic effect** across all levels of the subgrouping variable $j$.

- $\mathbf{v}_{i}$: A vector of covariates for which there is a strong
  *a priori* reason to believe their prognostic effect **varies** across
  the levels $k$ of the subgrouping variable $j$.

The linear predictor for patient $i$ in level $k$ of variable $j$ is
then:
$$g\left( \mu_{i,j,k} \right) = \beta_{1,j,k} + \beta_{2,j,k}z_{i} + {\mathbf{β}}_{3,j}\mathbf{u}_{i} + {\mathbf{β}}_{4,j,k}\mathbf{v}_{i}$$

- $\beta_{2,j,k}$: The predictive effect for level $k$ of variable $j$.
  This is the primary parameter of interest for shrinkage.
- ${\mathbf{β}}_{3,j}\mathbf{u}_{i}$: The effect of covariates assumed
  to be **constant** across levels. The coefficient vector
  ${\mathbf{β}}_{3,j}$ does **not** have a $k$ subscript.
- ${\mathbf{β}}_{4,j,k}\mathbf{v}_{i}$: The effect of covariates assumed
  to be **level-specific**. The coefficient vector
  ${\mathbf{β}}_{4,j,k}$**does** have a $k$ subscript, allowing its
  effect to be different for each subgroup level. This can be seen as an
  interaction between variable $x_{j}$ and $v_{i}$.

Note that for coefficients ${\mathbf{β}}_{3,j}$ and
${\mathbf{β}}_{4,j,k}$ we can use shrinkage or standard priors depending
in our prior beliefs.

#### 3.4.2 Example

Similarly to the example showed for the Global model. The example below
is for the **Age Group** variable, where
$k \in \{\text{<65},\text{≥65}\}$.

The outcome for a patient *i* in age subgroup *k* is:
$Y_{ik} \sim N\left( \mu_{ik},\sigma_{k}^{2} \right)$.

The linear predictor $\mu_{ik}$ is:

$$\mu_{ik} = \underset{\text{Subgroup Intercept}}{\underbrace{\beta_{1,k}}} + \underset{\text{Predictive Effect}}{\underbrace{\beta_{2,k} \cdot \text{Treatment}_{i}}} + \underset{\text{Constant Prognostic Effects}}{\underbrace{\delta_{1} \cdot \text{BMI}_{i} + \delta_{2} \cdot \text{Smoker}_{i}}} + \underset{\text{Varying Prognostic Effect}}{\underbrace{\delta_{3,k} \cdot \text{BaselineSBP}_{i}}}$$

- Here, $\beta_{1,k}$ (subgroup-specific intercept) and $\beta_{2,k}$
  (subgroup-specific predictive effect) are given hierarchical priors to
  borrow information across age groups. For example:
  $\beta_{2,k} \sim N\left( \mu_{2},\tau_{2}^{2} \right)$.
- The effect of Baseline SBP is also allowed to vary by age group
  ($\delta_{3,k}$), while BMI and Smoking effects are assumed constant

### 3.5 Endpoint-Specific Models

The package supports four primary endpoint types by specifying the
appropriate likelihood and link function for
$g({\mathbf{μ}}) = {\mathbf{η}}$.

| Endpoint          | Distribution / Model                                                       | Link Function                                        | Notes                                                                                                  |
|:------------------|:---------------------------------------------------------------------------|:-----------------------------------------------------|:-------------------------------------------------------------------------------------------------------|
| **Continuous**    | Normal: $y_{i} \sim N\left( \mu_{i},\sigma^{2} \right)$                    | Identity: $\mu_{i} = \eta_{i}$                       | $\sigma$ can be stratified (e.g., `sigma ~ 0 + clinic_site`).                                          |
| **Binary**        | Bernoulli: $y_{i} \sim \text{Bernoulli}\left( p_{i} \right)$               | Logit: $\text{logit}\left( p_{i} \right) = \eta_{i}$ |                                                                                                        |
| **Count**         | Negative Binomial: $y_{i} \sim \text{NegBin}\left( \mu_{i},\phi \right)$   | Log: $\log\left( \mu_{i} \right) = \eta_{i}$         | $\phi$ (overdispersion) can be stratified. Supports [`offset()`](https://rdrr.io/r/stats/offset.html). |
| **Time-to-event** | Cox Proportional Hazards: $h_{i}(t) = h_{0}(t)\exp\left( \eta_{i} \right)$ | Log (on the hazard)                                  | Intercept is absorbed by the baseline hazard $h_{0}(t)$. $h_{0}(t)$ can be stratified.                 |

### 3.6 Prior Specifications

#### 3.6.1 Weakly Informative Priors (Unshrunk Terms)

For all non-penalized coefficients
(${\mathbf{β}}_{\mathbf{1},\mathbf{n}},{\mathbf{β}}_{\mathbf{2},\mathbf{n}},{\mathbf{β}}_{\mathbf{3}}$),
we use weakly informative priors to aid computational stability without
strongly influencing the posterior.

- **Default:** `normal(0, 10)`. This is a common, weakly regularizing
  prior on the linear predictor scale. Users can specify their own,
  e.g., `student_t(3, 0, 2.5)`.

#### 3.6.2 Regularized Horseshoe Prior (Shrunk Terms)

This is the default shrinkage prior in `bonsaiforest2`, recommended for
its excellent adaptive shrinkage properties (Wolbers et al. 2025;
Piironen and Vehtari 2017).

- **Concept:** A global-local prior. It has an infinitely tall spike at
  zero (to aggressively shrink noise) and heavy tails (to leave true,
  large signals unshrunk) .
- **Hierarchical Specification** :

$$\begin{aligned}
\beta_{2,l} & {\sim N\left( 0,\tau^{2}{\widetilde{\lambda}}_{l}^{2} \right)} \\
{\widetilde{\lambda}}_{l}^{2} & {= \frac{c^{2}\lambda_{l}^{2}}{c^{2} + \tau^{2}\lambda_{l}^{2}}} \\
\lambda_{l} & {\sim C^{+}(0,1)\quad\left( \text{Local shrinkage} \right)} \\
\tau & {\sim C^{+}\left( 0,\tau_{0} \right)\quad\left( \text{Global shrinkage} \right)} \\
c^{2} & {\sim \text{Inverse-Gamma}\left( \nu/2,\nu s^{2}/2 \right)}
\end{aligned}$$

- **Hyperparameter Justification**: The package supports two ways to set
  the crucial global scale $\tau_{0}$:
  1.  **Fixed Default (`scale_global`):** `horseshoe(scale_global = 1)`.
      This is the package default, a robust, general-purpose choice
      (Wolbers et al. 2025).
  2.  **Elicited Prior (`par_ratio`):** Sets $\tau_{0}$ based on the
      *prior* belief about the number of effective (non-zero) subgroups,
      $p_{eff}$, out of the total $L$(Bornkamp 2025; Piironen and
      Vehtari 2017). This is specified via
      $\text{par\_ratio} = \frac{p_{eff}}{L - p_{eff}}$. This scales the
      prior based on sample size ($N$) and the expected sparsity.

#### 3.6.3 R2D2 Prior (Shrunk Terms)

- **Concept:** This prior is uniquely derived by first placing a prior
  on the model’s coefficient of determination, $R^{2}$, which quantifies
  the proportion of variance explained by the predictors. This is
  arguably more intuitive than specifying priors directly on
  coefficients.

- **Hierarchical Specification (Marginal Version):** $$\begin{aligned}
  \beta_{2,l} & {\sim DE\left( \sigma\sqrt{\phi_{l}\omega/2} \right)\quad\left( \text{Double-Exponential/Laplace kernel} \right)} \\
  {\mathbf{ϕ}} & {= \left( \phi_{1},\ldots,\phi_{L} \right) \sim \text{Dirichlet}\left( a_{\pi},\ldots,a_{\pi} \right)\quad\left( \text{Local shrinkage} \right)} \\
  \omega & {\sim \text{Beta-Prime}(a,b)\quad\left( \text{Global shrinkage} \right)}
  \end{aligned}$$ This is equivalent to assuming that the proportion of
  variance explained by the interaction terms follows
  $R^{2} = \frac{\omega}{1 + \omega} \sim \text{Beta}(a,b)$.

- **Justification of Hyperparameters:**

  - Zhang et al. (2022) provide clear guidance. A fully automatic
    approach is to set $b = 0.5$ to achieve Cauchy-like heavy tails.

  - The concentration parameter $a_{\pi}$ controls sparsity. A smaller
    $a_{\pi}$ (e.g., 0.2) implies higher shrinkage, concentrating the
    prior variance on fewer coefficients, while a larger $a_{\pi}$
    (e.g., 0.5) spreads the variance more evenly, implying lower
    shrinkage.

  - **Default Recommendation:** Set $b = 0.5$ and offer options for the
    user based on desired shrinkage strength:

    - **High Shrinkage:** $a_{\pi} = 0.2$

    - **Low Shrinkage:** $a_{\pi} = 0.5$

#### 3.6.4. Normal Hierarchical Prior (for One-Variable-at-a-Time Model)

This is the standard hierarchical model for subgroup analysis, which
assumes that treatment effects for levels within a subgrouping variable
are drawn from a common normal distribution. It’s less complex than the
Horseshoe or R2D2 priors but provides effective shrinkage, especially
when the number of subgroup levels is small.

- **Hierarchical Specification:** For a given subgrouping variable $j$
  with levels $k = 1,\ldots,K_{j}$:

$$\begin{aligned}
\beta_{2,j,k} & {\sim \mathcal{N}\left( \mu_{2,j},\tau_{2,j}^{2} \right)\quad\left( \text{Level-specific treatment effects} \right)} \\
\mu_{2,j} & {\sim \mathcal{N}\left( 0,s^{2} \right)\quad\left( \text{Common mean effect} \right)} \\
\tau_{2,j} & {\sim \text{Half-Normal}(a)\quad\left( \text{Between-subgroup heterogeneity} \right)}
\end{aligned}$$

- **Justification of Hyperparameters:** The key is setting the scale
  parameter a for the Half-Normal prior on the between-subgroup standard
  deviation, $\tau_{2,j}$, as this controls the degree of shrinkage. As
  recommended by Bornkamp (2025), this can be linked to the planned
  treatment effect ( $\delta_{plan}$) from the trial protocol, making
  the prior choice interpretable and consistent across different
  endpoints. The choice of a implies a prior on the plausible difference
  in treatment effects between any two subgroups.

### 3.7 Estimation (MCMC)

The joint posterior distribution of all parameters is complex and has no
closed-form solution. We use Markov Chain Monte Carlo (MCMC) methods to
draw samples from this distribution. Specifically, the package uses the
**No-U-Turn Sampler (NUTS)**, an efficient variant of Hamiltonian Monte
Carlo (HMC), as implemented in **Stan** via the `brms` package. The
output is a set of posterior samples (e.g., 4000 draws) representing the
joint posterior distribution.

### 3.8 Standardization for Marginal Effects (G-computation)

To obtain interpretable *marginal* effects for each subgroup, the
package implements a standardization (G-computation) procedure. This
process is repeated for *each* MCMC draw to generate a full posterior
distribution of the marginal effect.

**For a single MCMC draw, the step-by-step process is:**

1.  **Select Parameters:** Take one draw from the joint posterior
    distribution for all model parameters (all $\mathbf{β}$’s, $\sigma$,
    $\phi$, etc.).

2.  **Create Counterfactuals:** For *every patient* $i$ in the original
    dataset, create two scenarios:

    - Scenario A: Patient $i$’s covariates, with treatment $s_{i} = 0$
      (Control).
    - Scenario B: Patient $i$’s covariates, with treatment $s_{i} = 1$
      (Treatment).

3.  **Predict Outcomes:** Use the model formula and the parameters from
    Step 1 to predict the outcome for *every patient* under both
    Scenario A (${\widehat{\mu}}_{i,0}$) and Scenario B
    (${\widehat{\mu}}_{i,1}$).

4.  **Average within Subgroups:** For a subgroup of interest $l$ (e.g.,
    “Female”):

    - Calculate the average predicted outcome under control:
      ${\widehat{\mu}}_{l,0} = \text{mean}\left( {\widehat{\mu}}_{i,0} \right)$
      for all $i$ where $x_{il} = 1$.
    - Calculate the average predicted outcome under treatment:
      ${\widehat{\mu}}_{l,1} = \text{mean}\left( {\widehat{\mu}}_{i,1} \right)$
      for all $i$ where $x_{il} = 1$.

5.  **Calculate Effect Measure:** Compute the marginal effect for
    subgroup $l$ from these averaged predictions. This depends on the
    endpoint type:

    - **Continuous:** Mean Difference =
      ${\widehat{\mu}}_{l,1} - {\widehat{\mu}}_{l,0}$.
    - **Binary:** Odds Ratio =
      $\frac{{\widehat{\mu}}_{l,1}/\left( 1 - {\widehat{\mu}}_{l,1} \right)}{{\widehat{\mu}}_{l,0}/\left( 1 - {\widehat{\mu}}_{l,0} \right)}$
      (where $\widehat{\mu}$ is the predicted probability).
    - **Count:** Rate Ratio =
      ${\widehat{\mu}}_{l,1}/{\widehat{\mu}}_{l,0}$ (where
      $\widehat{\mu}$ is the predicted rate).
    - **Time-to-Event:** Average Hazard Ratio (AHR). This is the most
      complex, as it requires predicting the full survival curve
      $S_{i}(t)$ for each patient. The marginal survival curves
      ${\widehat{S}}_{l,0}(t)$ and ${\widehat{S}}_{l,1}(t)$ are
      calculated, and the AHR is computed from them, as the marginal
      curves may not be proportional.

6.  **Repeat:** This process (Steps 1-5) is repeated for all MCMC draws,
    yielding a full posterior distribution (e.g., 4000 values) for the
    marginal effect of each subgroup. The final point estimate (median)
    and 95% credible interval are taken from this distribution.

## 4. Mapping of Statistical Methods to `bonsaiforest2` Functions

This section connects the statistical methodology (Section 3) to the
core package functions.

- **[`run_brms_analysis()`](https://openpharma.github.io/bonsaiforest2/reference/run_brms_analysis.md)**

  - **Maps to:** Sections 3.2 (Global Model), 3.3 (Endpoints), 3.4
    (Priors), 3.5 (Estimation).
  - **Action:** This function is the primary model-fitting engine.
    - It constructs the design matrices
      $\mathbf{X}_{\mathbf{n}},\mathbf{X}_{\mathbf{p}},\mathbf{Z}_{\mathbf{n}},\mathbf{Z}_{\mathbf{p}},\mathbf{U}$
      based on the formula string arguments (e.g.,
      `unshrunk_prognostic_formula_str`,
      `shrunk_predictive_formula_str`).
    - It assigns the specified priors (e.g.,
      `predictive_effect_priors = list(shrunk = "horseshoe(par_ratio = 0.1)")`).
    - It builds and passes the full `brmsformula` and data to
      [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html)
      to run the Stan MCMC sampler.
    - It handles stratification via the `stratification_formula_str`
      argument, which adds terms like `bhaz(country)` or
      `sigma ~ country` to the `brmsformula`.

- **[`summary_subgroup_effects()`](https://openpharma.github.io/bonsaiforest2/reference/summary_subgroup_effects.md)**

  - **Maps to:** Section 3.6 (Standardization / G-computation).
  - **Action:** This function implements the full G-computation
    procedure.
    - It takes the fitted `brms_fit` object (containing posterior
      samples) and the `original_data`.
    - It automatically detects the subgroups specified in the model (or
      uses the `subgroup_vars` argument).
    - It iterates through each posterior sample, generates
      counterfactual predictions, averages within subgroups, and
      calculates the marginal effect measure (e.g., AHR, OR) based on
      the `response_type`.
    - It returns a `subgroup_summary` object containing the posterior
      distributions of the marginal effects.

- **[`plot()`](https://rdrr.io/r/graphics/plot.default.html)**

  - **Maps to:** Final visualization of Section 3.6 results.
  - **Action:** This is a method for `subgroup_summary` objects. It
    takes the summarized posterior distributions (median and 95% CrI)
    for the overall and subgroup effects and generates a forest plot.

## 5. References

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

------------------------------------------------------------------------

1.  Note that this is different from what they do in the Wang et
    al. (2024) where they give shrinkage priors to all the coefficients.
    This can also be done in our library.
