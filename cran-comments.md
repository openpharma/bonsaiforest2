## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

* Suggested package 'cmdstanr' is not on CRAN. It is available from
  https://stan-dev.r-universe.dev, declared via `Additional_repositories`.
  It is optional: models are fitted with 'rstan' when 'cmdstanr' is not
  installed. Tests and vignettes only use it when available, and the
  one example that requests it is wrapped in `\dontrun{}`.

## Test environments

* local: macOS Tahoe 26.5.2 (aarch64), R 4.6.1
* win-builder: R-devel (2026-11-02 r90827 ucrt)
* win-builder: R-release 4.6.2 (2026-10-29 ucrt)
