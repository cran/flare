# flare

`flare` provides a family of sparse regression and sparse graphical model estimators,
including Dantzig Selector, LAD Lasso, SQRT Lasso, Lq Lasso, TIGER, and CLIME.

## Development Notes

- The CI workflow runs:
  - `R CMD check` on Linux/macOS/Windows (R release)
  - `R CMD check` on Linux (R devel)
  - Test coverage generation and optional Codecov upload
