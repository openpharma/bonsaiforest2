# Statistical Specifications

**A comprehensive description of the statistical methodology is provided
in Wolbers et al. ([2026](#ref-wolbersUnifiedShrinkage)). In this
vignette, we only provide an overview.**

## 1 Introduction

The `bonsaiforest2` package implements Bayesian shrinkage models for
treatment effect estimation in subgroups or randomized clinical trials.
It supports two approaches to Bayesian shrinkage: **One-way shrinkage
models** fit a Bayesian hierarchical model for one subgrouping variable
at a time (Wang et al. ([2024](#ref-wang2024bayesian))), typically
requiring multiple model fits to obtain effects for all subgroups of
interest. **Global shrinkage models** fit a single model including
treatment interactions for all subgroups jointly (Wolbers et al.
([2025](#ref-wolbers2025using))).

The estimation process consists of two steps:

1.  **Model Fitting:** Fit a Bayesian model to the trial data, i.e. a
    one-way or a global model
2.  **Standardization:** Derive standardized treatment effect estimates
    for subgroups via G-computation (marginalizing over the covariate
    distribution within the subgroup)

## 2 Model Specification and Fitting Step

### 2.1 Notation

We denote the subjects included in the trial by the subscript \\i\\ (\\i
= 1,\ldots , N\\) and their treatment assignment by \\z_i\\ taking
values of 0 for control and 1 for intervention. For the **one-way
shrinkage model**, assume that we are interested in treatment effects in
subgroups defined by the levels of a single categorical subgrouping
variable \\x\\ with \\K\\ levels. For **global shrinkage models**,
assume that subgroups are deﬁned by the levels of \\p\\ categorical
subgrouping variables \\x_j\\ \\(j = 1,\ldots , p)\\, respectively. Each
variable \\x_j\\ has \\l_j\\ levels, resulting in a total of \\K =
\sum\_{j=1}^p l_j\\ subgroups defined by a single variable at a time.
For both models, denote indicators for each of the \\K\\ subgroups by
\\s\_{ik}\\ (i.e., \\s\_{ik}=1\\ if subject \\i\\ belongs to the
subgroup \\k\\, and \\s\_{ik}=0\\ otherwise, with \\k=1,\ldots, K\\). In
addition, the model might be adjusted for additional baseline covariates
denoted by \\u\_{il}\\ (\\l=1,\ldots,L\\).

| Feature                 | One-Way shrinkage Model                                       | Global shrinkage Model                                                             |
|:------------------------|:--------------------------------------------------------------|:-----------------------------------------------------------------------------------|
| **Scope**               | Fit one model per subgrouping variables                       | Fit one model including all the subgrouping variables                              |
| **Subgroups (\\K\\)**   | \\K\\ levels of a single subgrouping variable (disjoint)      | \\K = \sum\_{j=1}^p l_j\\ subgroups from \\p\\ subgrouping variables (overlapping) |
| **Subgroup Indicators** | Only a single indicator \\s\_{ik}\\ is equal to 1 per subject | \\p\\ indicators \\s\_{ik}\\ are equal to 1 per subject                            |

### 2.2 Statistical model

The outcome model is a Bayesian regression model with a linear predictor
structure. `bonsaiforest2` supports multiple endpoint types and
treatment effect measures:

- **Continuous:** Linear regression (Mean differences)
- **Binary:** Logistic regression (Odds Ratios, OR)
- **Count:** Negative Binomial regression (Rate Ratios, RR) — *supports
  offset terms*
- **Time-to-Event:** Cox regression (Average Hazard Ratios, HR) — *uses
  non-negative M-splines for the baseline hazard*

Regardless of the endpoint type, the linear predictor \\LP_i\\ for
subject \\i\\ is structured as follows:

\\LP_i = \alpha_0 + \beta\_{0}z_i + \underbrace{\alpha_1
s\_{i1}+\ldots + \alpha_K s\_{iK} + \alpha\_{K+1} u\_{i1} +\ldots +
\alpha\_{K+L}u\_{iL} }\_{\substack{\text{prognostic effects including}
\\ \text{main subgroup effects}}} + \underbrace{\beta_1 s\_{i1} z_i
+\ldots + \beta_K s\_{iK}z_i}\_{\substack{\text{predictive effects:
subgroup-by-} \\ \text{treatment interactions}}}\\

*Note: For Cox regression models, the intercept \\\alpha_0\\ is
excluded.*

### 2.3 Parameter grouping

To provide model flexibility, `bonsaiforest2` allows the user to
categorize regression coefficients of the linear predictor (excluding
the intercept \\\alpha_0\\ and main treatment effect \\\beta_0\\) into
three distinct groups.

| Parameter Group            | **Unshrunken Terms**                                                                                                                                                                                                                                                                                                   | **Shrunken Prognostic Terms**                                                                                                                                                                                                                     | **Shrunken Predictive Terms**                                                                                                                                                                                                                                                     |
|:---------------------------|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Purpose**                | Typically consists of all prognostic effects. If there is a strong a priori rationale for treatment effect heterogeneity for some of the subgroups in a global model, e.g., for those defined by biomarkers directly linked to a drug target, then they might also be left unshrunken to reflect this prior knowledge. | Prognostic terms with shrinkage priors. Often, this group is empty as all prognostic terms are left unshrunken. Prognostic terms might be shrunken if their number is large (relative to the sample size) and regularization is deemed necessary. | For one-way shrinkage models, this group includes the subgroup-by-treatment interactions for the selected subgrouping variable only. For global models, it encompasses all subgroup-by-treatment interactions (with the possible exception of those grouped as unshrunken terms). |
| **Desing matrix encoding** | Dummy encoding (baseline levels omitted to avoid overparameterization).                                                                                                                                                                                                                                                | One-hot encoding (all levels included symmetrically)                                                                                                                                                                                              | One-hot encoding (all levels included symmetrically)                                                                                                                                                                                                                              |
| **Prior distribution**     | Non-informative / Flat priors                                                                                                                                                                                                                                                                                          | Shrinkage priors (exchangeable within this group)                                                                                                                                                                                                 | Shrinkage priors (exchangeable within this group)                                                                                                                                                                                                                                 |

### 2.4 Priors distributions for shrunken predictive terms

Our implementation is flexible in terms of prior choices. A common
choice for one-way shrinkage models (where the number of subgroups is
typically low, e.g. 2-5) are normal priors. A common choice for global
shrinkage models (where the number of subgroups is larger, e.g. 10-30)
are regularized horseshoe priors.

#### 2.4.1 Normal priors for one-way shrinkage models

A normal prior distribution with a half-normal (HN) hyperprior assigned
to the standard deviation is the typical choice:

\\\beta_1, \ldots, \beta_K \sim N(0,\tau^2) \mbox{ with } \tau \sim
HN(\phi)\\

The heterogeneity parameter \\\phi\\ dictates the degree of shrinkage. A
possible **default choice** is \\\phi = \|\delta\_{plan}\|\\ where
\\\delta\_{plan}\\ is the target effect size from the trial protocol.
This choice implies an a priori assumption that pairwise differences in
treatment effects between subgroups (i.e. \\\|\beta_i-\beta_j\|\\) lie
between \\0.02\cdot\|\delta\_{plan}\|\\ and \\3.09\cdot
\|\delta\_{plan}\|\\ with 90% probability. This range seems sufficiently
wide to remain conservative, even in settings characterized by
substantial heterogeneity.

#### 2.4.2 Regularized horseshoe priors for global shrinkage models

The larger number of subgroups encourages the use of global-local
shrinkage priors. Global-local shrinkage priors are continuous versions
of spike-and-slab priors, i.e., mixture priors of a point mass at zero
(the spike) and a continuous distribution which describes the non-zero
coefficients (the slab). A common choice is the regularized horseshoe
prior which has the following form:

\\ \beta_k\|\lambda_k, \tau, c \sim
N(0,\tau^2\tilde{\lambda}\_k^2)\mbox{ with } \tilde{\lambda}\_k^2 = c^2
\lambda_k^2 / (c^2 + \tau^2 \lambda_k^2)\\

with hyper-priors: \\\tau\sim C^+(0,\tau_0^2), \\ \lambda_k \sim
C^+(0,1)\\ (k=1,\ldots,K), \mbox{ and } c^2 \sim
\textrm{Inv-Gamma}(\nu/2,\nu s^2 /2)\\ The global shrinkage parameter
\\\tau\\ shrinks all parameters towards zero while the heavy-tailed
half-Cauchy prior for \\\lambda_k\\ allows some parameters to escape
this heavy shrinkage. The slab-component for the regularized horseshoe
which governs shrinkage for larger coefficients is a
\\Student-t\_{\nu}(0,s^2)\\ distribution (Piironen and Vehtari
([2017](#ref-piironen2017))).

As described in Wolbers et al. ([2026](#ref-wolbersUnifiedShrinkage)), a
possible **default choice** is \\\tau_0 = \|\delta\_{plan}\|\\
(`scale_global`), \\s = 2\sigma\_{plan}\\ for continuous and \\s=2\\ for
other endpoints (`scale_slab`), and \\\nu = 4\\ (`df_slab`).

### 2.5 Model fitting

In our implementation, we specify the outcome model using the `R`
package `brms`, which relies on Stan as the backend for NUTS or
Hamiltonian Monte Carlo sampling (Bürkner ([2017](#ref-bruckner2017)),
Carpenter et al. ([2017](#ref-carpenter2017stan))).

## 3 Standardization Step (G-Computation)

The coefficients \\\beta_i\\ from the regression model cannot be
interpreted as standardized treatment effects in subgroups.
`bonsaiforest2` uses standardization to convert them into clinically
interpretable treatment effects.

The standardization process follows three steps:

1.  **Prediction:** The fitted regression model is used to predicts
    potential outcomes for every subject \\i\\ twice: once assuming
    hypothetical assignment to control (\\z_i=0\\) and once to
    intervention (\\z_i=1\\). For time-to-event data, full survival
    curves \\\hat{S}\_{i,C}(t)\\ and \\\hat{S}\_{i,I}(t)\\ are
    predicted.
2.  **Marginalization:** Predicted potential outcomes (or survival
    curves) are averaged across all subjects *within* a specific
    subgroup, resulting in a marginal prediction for the subgroup under
    both treatment conditions.
3.  **Estimation:** The subgroup treatment effect is calculated by
    contrasting the marginal intervention prediction against the
    marginal control prediction (e.g., yielding a marginal mean
    difference, OR, RR, or average HR).

Standardized treatment effect estimates in subgroups are summarized
using the median of the posterior draws, accompanied by a 95% credible
interval (bounded by the 2.5% and 97.5% quantiles).

## 4 Comparison of One-way vs Global Shrinkage Models

To us, global shrinkage models are a more principled approach than
one-way models to obtain treatment effect estimates across a
well-defined family of multiple subgrouping variables.

In terms of actual performance, the following high-level conclusions may
be derived from the simulation studies reported in Wolbers et al.
([2026](#ref-wolbersUnifiedShrinkage)):

- Both one-way and global shrinkage models typically have a
  substantially lower overall root mean squared error than the standard
  estimator without shrinkage. They are also more efficient in
  identifying a non-efficacious subgroup.
- Global shrinkage models tend to have lower mean-squared error and are
  less dependent on the choice of prior parameters.
- One-way shrinkage models tend to have slightly less bias and better
  frequentist coverage of associated credible intervals.
- Covariate-adjustment for known prognostic factors can substantially
  improve the precision of treatment effect estimates in subgroups.
- For both models, hyperprior parameter choices anchored in trial
  assumptions about the anticipated size of the overall treatment effect
  (such as the possible default choices described above) performed well.

## References

Bürkner, Paul-Christian. 2017. “brms: An R Package for Bayesian
Multilevel Models Using Stan.” *Journal of Statistical Software* 80 (1):
1–28.

Carpenter, Bob, Andrew Gelman, Matthew D Hoffman, Daniel Lee, Ben
Goodrich, Michael Betancourt, Marcus Brubaker, Jiqiang Guo, Peter Li,
and Allen Riddell. 2017. “Stan: A Probabilistic Programming Language.”
*Journal of Statistical Software* 76: 1–32.

Piironen, Juho, and Aki Vehtari. 2017. “Sparsity Information and
Regularization in the Horseshoe and Other Shrinkage Priors.” *Electronic
Journal of Statistics* 11 (2): 5018–51.

Wang, Yun, Wenda Tu, William Koh, James Travis, Robert Abugov, Kiya
Hamilton, Mengjie Zheng, Roberto Crackel, Pablo Bonangelino, and Mark
Rothmann. 2024. “Bayesian hierarchical models for subgroup analysis.”
*Pharmaceutical Statistics* 23: 1065–83.

Wolbers, Marcel, Miriam Pedrera Gómez, Alex Ocampo, and Isaac
Gravestock. 2026. “Unified Implementation and Comparison of Bayesian
Shrinkage Methods for Treatment Effect Estimation in Subgroups.” *arXiv
Preprint arXiv:2603.21967*. <https://arxiv.org/abs/2603.21967>.

Wolbers, Marcel, Mar Vázquez Rabuñal, Ke Li, Kaspar Rufibach, and Daniel
Sabanés Bové. 2025. “Using shrinkage methods to estimate treatment
effects in overlapping subgroups in randomized clinical trials with a
time-to-event endpoint.” *Statistical Methods in Medical Research*,
1–12.
