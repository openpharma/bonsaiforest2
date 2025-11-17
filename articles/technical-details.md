# Technical Details of bonsai forest Models

## Introduction

This vignette details the statistical framework, notation, models, and
prior specifications implemented within the `bonsaiforest2` package. It
is intended for users seeking a deeper understanding of the underlying
methodology used for Bayesian subgroup analysis.

For a practical guide on how to use the package functions, please see
the main user guide vignette:
[`vignette("using-bonsaiforest", package = "bonsaiforest2")`](https://openpharma.github.io/bonsaiforest2/articles/using-bonsaiforest.md).

## General Model Framework

### Notation

We establish the following notation:

- $i = 1,\ldots,N$: Index for individual patients.
- $y_{i}$: The outcome for patient $i$.
- $s_{i}$: The binary treatment indicator for patient $i$ (e.g.,
  $s_{i} = 1$ for treatment, $s_{i} = 0$ for control).
- $j = 1,\ldots,J$: Index for the subgrouping *variable* (e.g., age
  group, region).
- $k = 1,\ldots,K_{j}$: Index for the subgroup *level* within variable
  $j$.
- $l = 1,\ldots,L$: Index for an overall subgroup level across all
  variables ($L = \sum_{j = 1}^{J}K_{j}$).
- $x_{il}$: Indicator variable (0/1) for patient $i$ belonging to
  overall subgroup level $l$.
- $z_{il} = s_{i} \cdot x_{il}$: The interaction between treatment and
  subgroup level $l$ for patient $i$.
- $u_{iv}$: Value of the $v$-th additional non-subgroup covariate for
  patient $i$.

### Prognostic vs. Predictive Effects

In clinical trials, baseline characteristics can relate to outcomes in
two ways:

- **Prognostic Effect**: Association with the outcome, irrespective of
  treatment (a **main effect**).
- **Predictive Effect**: Modifies the treatment effect (a
  **treatment-by-variable interaction**).

`bonsaiforest2` explicitly models both types of effects.

### Design Matrix Considerations

The package constructs design matrices differently for prognostic and
predictive effects based on modeling goals:

- **Prognostic Effects (`X` matrix):** Uses standard **dummy coding**
  (dropping one reference level per categorical variable). This avoids
  collinearity with the intercept and yields interpretable main effect
  coefficients relative to a baseline group.
- **Predictive Effects (`Z` matrix):** Uses **one-hot encoding**
  (creating a binary indicator for *every* level of each categorical
  variable, *without* dropping a reference). This treats all subgroup
  levels symmetrically under the shrinkage prior, reflecting an
  assumption of exchangeability for the interaction effects.

### General Global Model

The package primarily implements a global model structure (Wolbers et
al. 2025) where all effects are estimated simultaneously. The linear
predictor, $\mathbf{η}$, for all $N$ patients is:

$$g({\mathbf{μ}}) = {\mathbf{η}} = \mathbf{X}_{\mathbf{p}}{\mathbf{β}}_{\mathbf{1},\mathbf{p}} + \mathbf{X}_{\mathbf{n}}{\mathbf{β}}_{\mathbf{1},\mathbf{n}} + \mathbf{Z}_{\mathbf{p}}{\mathbf{β}}_{\mathbf{2},\mathbf{p}} + \mathbf{Z}_{\mathbf{n}}{\mathbf{β}}_{\mathbf{2},\mathbf{n}} + \mathbf{U}{\mathbf{β}}_{\mathbf{3}}$$
Where:

- $g( \cdot )$ is the link function.
- $\mathbf{μ}$ is the vector of expected outcomes.
- $\mathbf{X}_{\mathbf{n}}$ contains columns for **unshrunk prognostic**
  effects (intercept, main treatment effect, specified main effects)
  using dummy coding. ${\mathbf{β}}_{\mathbf{1},\mathbf{n}}$ are their
  coefficients.
- $\mathbf{X}_{\mathbf{p}}$ contains columns for **shrunk prognostic**
  effects (specified main effects) using dummy coding.
  ${\mathbf{β}}_{\mathbf{1},\mathbf{p}}$ are their coefficients,
  assigned a shrinkage prior.
- $\mathbf{Z}_{\mathbf{n}}$ contains columns for **unshrunk predictive**
  effects (treatment interactions) using one-hot encoding.
  ${\mathbf{β}}_{\mathbf{2},\mathbf{n}}$ are their coefficients.
- $\mathbf{Z}_{\mathbf{p}}$ contains columns for **shrunk predictive**
  effects (treatment interactions) using one-hot encoding.
  ${\mathbf{β}}_{\mathbf{2},\mathbf{p}}$ are their coefficients,
  assigned a shrinkage prior.
- $\mathbf{U}$ contains columns for additional **unshrunk prognostic**
  covariates. ${\mathbf{β}}_{\mathbf{3}}$ are their coefficients.

This structure allows flexible assignment of terms to shrunk or unshrunk
categories:

|                       | **Shrinkage Prior**                                           | **Weakly Informative Prior**                                                                          |
|:----------------------|:--------------------------------------------------------------|:------------------------------------------------------------------------------------------------------|
| **Prognostic Effect** | $\mathbf{X}_{\mathbf{p}}{\mathbf{β}}_{\mathbf{1},\mathbf{p}}$ | $\mathbf{X}_{\mathbf{n}}{\mathbf{β}}_{\mathbf{1},\mathbf{n}}$ + $\mathbf{U}{\mathbf{β}}_{\mathbf{3}}$ |
| **Predictive Effect** | $\mathbf{Z}_{\mathbf{p}}{\mathbf{β}}_{\mathbf{2},\mathbf{p}}$ | $\mathbf{Z}_{\mathbf{n}}{\mathbf{β}}_{\mathbf{2},\mathbf{n}}$                                         |

The linear predictor for a single patient $i$ is:

$$\eta_{i} = \underset{\text{Prognostic Component}}{\underbrace{\beta_{0} + \beta_{\text{treat}}s_{i} + \sum\beta^{\text{prog, shrunk}}x_{i} + \sum\beta^{\text{prog, unshrunk}}x_{i} + \sum\beta^{\text{extra}}u_{i}}} + \underset{\text{Predictive Component}}{\underbrace{\sum\beta^{\text{pred, shrunk}}z_{i} + \sum\beta^{\text{pred, unshrunk}}z_{i}}}$$*(Note:
The indices and notation are simplified here for clarity compared to the
full matrix form).*

## Endpoint-Specific Models

The package supports four endpoint types by selecting the appropriate
likelihood and link function within `brms`.

| Endpoint          | Distribution / Model                                                       | Link Function                                                                                     | Notes                                                                                                         |
|:------------------|:---------------------------------------------------------------------------|:--------------------------------------------------------------------------------------------------|:--------------------------------------------------------------------------------------------------------------|
| **Continuous**    | Normal: $y_{i} \sim N\left( \mu_{i},\sigma^{2} \right)$                    | Identity: $\mu_{i} = \eta_{i}$                                                                    | Residual standard deviation $\sigma$ can be stratified.                                                       |
| **Binary**        | Bernoulli: $y_{i} \sim \text{Bernoulli}\left( p_{i} \right)$               | Logit: $\text{logit}\left( p_{i} \right) = \log\left( \frac{p_{i}}{1 - p_{i}} \right) = \eta_{i}$ |                                                                                                               |
| **Count**         | Negative Binomial: $y_{i} \sim \text{NegBin}\left( \mu_{i},\phi \right)$   | Log: $\log\left( \mu_{i} \right) = \eta_{i}$                                                      | Overdispersion $\phi$ can be stratified.                                                                      |
| **Time-to-event** | Cox Proportional Hazards: $h_{i}(t) = h_{0}(t)\exp\left( \eta_{i} \right)$ | Log (on the hazard)                                                                               | Intercept absorbed by baseline hazard $h_{0}(t)$, modeled via `bhaz()` splines. $h_{0}(t)$ can be stratified. |

------------------------------------------------------------------------

## Prior Specifications

Appropriate prior selection is crucial in Bayesian analysis, especially
when using shrinkage.

### Weakly Informative Priors

For parameters *not* subject to shrinkage (intercept, main treatment
effect, covariates in $\mathbf{X}_{\mathbf{n}}$,
$\mathbf{Z}_{\mathbf{n}}$, $\mathbf{U}$), **weakly informative priors**
are used. These provide minimal regularization, primarily to aid
computation and prevent unreasonable parameter estimates.

- **Default:** $N\left( 0,s^{2} \right)$, where $s$ is large but
  sensible (e.g., $s = 10$ on the linear predictor scale). The
  `bonsaiforest2` package uses `normal(0, 10)` as the default. Users can
  override this.

### Shrinkage Priors for Penalized Terms (${\mathbf{β}}_{\mathbf{1},\mathbf{p}},{\mathbf{β}}_{\mathbf{2},\mathbf{p}}$)

These priors pull coefficients towards a common center (usually zero),
effectively reducing noise and performing model selection.
`bonsaiforest2` facilitates the use of several state-of-the-art
shrinkage priors via `brms`.

*Marginal density comparison of various shrinkage priors. Reference:
Zhang et al. (2022).*

#### Regularized Horseshoe Prior

- **Concept:** A highly effective global-local prior with a sharp peak
  at zero (strong shrinkage for null effects) and heavy tails (minimal
  shrinkage for large, true effects) (Wolbers et al. 2025). The
  “regularized” version adds slight penalization for very large effects.

- **Specification:** $$\begin{aligned}
  \beta_{l} & {\sim N\left( 0,\tau^{2}{\widetilde{\lambda}}_{l}^{2} \right)} \\
  {\widetilde{\lambda}}_{l}^{2} & {= \frac{c^{2}\lambda_{l}^{2}}{c^{2} + \tau^{2}\lambda_{l}^{2}}} \\
  \lambda_{l} & {\sim C^{+}(0,1)\quad\left( \text{Local shrinkage} \right)} \\
  \tau & {\sim C^{+}\left( 0,\tau_{0} \right)\quad\left( \text{Global shrinkage} \right)} \\
  c^{2} & {\sim \text{Inverse-Gamma}\left( \nu/2,\nu s^{2}/2 \right)\quad\left( \text{Regularization} \right)}
  \end{aligned}$$

- **Hyperparameter Justification:** The global scale $\tau_{0}$ is key.
  `brms` (and `bonsaiforest2`) allows specification via
  [`horseshoe()`](https://paulbuerkner.com/brms/reference/horseshoe.html)
  arguments:

  - `scale_global`: Directly sets $\tau_{0}$. A common default is
    `scale_global = 1`.
  - `par_ratio`: Sets $\tau_{0}$ based on the expected proportion of
    non-zero coefficients ($p_{eff}$) out of the total ($L$), using
    $\tau_{0} = \left( p_{eff}/\left( L - p_{eff} \right) \right)/\sqrt{N}$(Piironen
    and Vehtari 2017; Bornkamp 2025). This allows tailoring shrinkage
    based on prior belief about heterogeneity.

- **Package Default:** The default shrinkage prior in `bonsaiforest2` is
  `horseshoe(scale_global = 1)`, following Wolbers et al. (2025). Users
  can specify alternatives (e.g., using `par_ratio`).

#### R2D2 Prior

- **Concept:** Places a prior on the proportion of variance explained
  ($R^{2}$) by the shrunk coefficients, then allocates this variance
  using a Dirichlet distribution, which controls sparsity (Zhang et al.
  2022; Bornkamp 2025).
- **Specification (Simplified):** Coefficients $\beta_{l}$ follow a
  distribution related to Double-Exponential (Laplace), where the
  individual variances depend on a global variance term $\omega$
  (related to $R^{2}$ via a Beta-Prime or Beta distribution) and local
  allocation parameters $\mathbf{ϕ}$ (from a Dirichlet distribution).

$$\begin{aligned}
\beta_{l} & {\sim \text{Kernel}\left( \sigma,\phi_{l},\omega \right)} \\
{\mathbf{ϕ}} & {\sim \text{Dirichlet}\left( a_{\pi},\ldots,a_{\pi} \right)} \\
\omega & {\sim \text{Beta-Prime}(a,b)\quad\left( \Leftrightarrow R^{2} = \frac{\omega}{1 + \omega} \sim \text{Beta}(a,b) \right)}
\end{aligned}$$ - **Hyperparameter Justification:** Following Zhang et
al. (2022): - Set $b = 0.5$ for heavy tails (recommended default). - The
Dirichlet concentration $a_{\pi}$ controls sparsity: smaller $a_{\pi}$
(e.g., 0.2) implies higher shrinkage (more sparsity); larger $a_{\pi}$
(e.g., 0.5 or 1.0) implies lower shrinkage. - **Package
Implementation:** Users can specify the R2D2 prior by providing the
string `"R2D2(mean_R2=..., prec_R2=..., cons_D=...)"` to the prior
arguments, where `cons_D` corresponds to $a_{\pi}$.

## Estimation and Inference Procedure

### Posterior Sampling: MCMC

Models are fitted using Hamiltonian Monte Carlo (HMC) via Stan, as
implemented in `brms`. This generates samples from the joint posterior
distribution of all model parameters.

### Deriving Subgroup Treatment Effects

For the global model, the treatment effect for subgroup level $l$ on the
linear predictor scale is
$\theta_{l} = \beta_{\text{treat}} + \beta_{2,l}$ (where $\beta_{2,l}$
is the corresponding interaction coefficient from either
${\mathbf{β}}_{\mathbf{2},\mathbf{p}}$ or
${\mathbf{β}}_{\mathbf{2},\mathbf{n}}$). This sum is computed for each
MCMC draw to obtain the posterior distribution of $\theta_{l}$.

### Standardization for Marginal Effects (G-computation)

Model coefficients typically represent *conditional* effects. To obtain
*marginal* effects (average treatment effect within a subgroup, averaged
over other covariates), the package uses G-computation implemented via
posterior predictions.

**Process per MCMC draw:**

1.  **Create Counterfactuals:** For each subject $i$, predict their
    outcome under treatment ($s_{i} = 1$) and control ($s_{i} = 0$),
    keeping all other covariates $u_{i},x_{i}$ as observed.
2.  **Predict Outcomes:** Use the parameter values from the current MCMC
    draw to calculate the predicted outcome (e.g., probability, rate,
    mean, survival curve) for each subject under both scenarios.
3.  **Average within Subgroups:** For a specific subgroup level $l$,
    calculate the average predicted outcome under treatment
    (${\widehat{\mu}}_{l,\text{treat}}$) and control
    (${\widehat{\mu}}_{l,\text{cont}}$) *only among subjects belonging
    to subgroup* $l$.
4.  **Calculate Effect Measure:** Compute the desired marginal effect
    based on the response type (e.g.,
    ${\widehat{\mu}}_{l,\text{treat}} - {\widehat{\mu}}_{l,\text{cont}}$
    for continuous; Odds Ratio for binary; Rate Ratio for count; Average
    Hazard Ratio for survival). This yields one draw from the posterior
    of the marginal effect for subgroup $l$.
5.  **Repeat:** Repeat for all MCMC draws to get the full posterior
    distribution for each subgroup’s marginal effect. Summaries (median,
    CrI) are then calculated from these distributions.

This standardization ensures that subgroup comparisons are made on a
comparable, marginal scale, accounting for differing covariate
distributions across subgroups.

## References

Bornkamp, Björn. 2025. “Benchmarking Bayesian subgroup shrinkage methods
on clinical data.” In *PSI Conference 2025*. London, UK.

Piironen, J., and A. Vehtari. 2017. “Sparsity information and
regularization in the horseshoe and other shrinkage priors.” *Electronic
Journal of Statistics* 11: 5018–51.

Wolbers, Marcel, Mar Vázquez Rabuñal, Ke Li, Kaspar Rufibach, and Daniel
Sabanés Bové. 2025. “Using shrinkage methods to estimate treatment
effects in overlapping subgroups in randomized clinical trials with a
time-to-event endpoint.” *Statistical Methods in Medical Research*,
1–12.

Zhang, Yan Dora, Brian P. Naughton, Howard D. Bondell, and Brian J.
Reich. 2022. “Bayesian regression using a prior on the model fit: The
R2-D2 shrinkage prior.” *Journal of the American Statistical
Association* 117 (538): 862–74.
