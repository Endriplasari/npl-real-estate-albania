# Change log

## 2.0.0 (revised manuscript, November 2026)

This version accompanies the revised manuscript submitted in response to the
reviewers' reports. The outputs of version 1.1.0, which accompanied the first
submission, are kept in `output/v1_submitted/` so that every change can be
traced.

### Corrections to the computations of version 1

| # | What was wrong in version 1 | Effect | Where it is fixed |
|---|---|---|---|
| 1 | Residual diagnostics were computed on `u.model.exog`, i.e. the static levels regression, not on the estimated error-correction model. | Table 8 of the first submission: identical rows for M1 and M1t and rejections of every test in every model. On the correct residuals RESET and CUSUM do not reject in any model. | `common.diagnostics` |
| 2 | The bounds F-statistic came from statsmodels' `UECM.bounds_test`, which drops fixed regressors, so the COVID dummy was excluded from the tested model. | F-statistics change (M1: 1.805 to 2.290). No decision changes. | `common.bounds_stats` |
| 3 | The asymptotic critical values reported were those for k + 1 regressors. | The "asymptotic 5%" column was too low (M1: 2.64 / 3.79 instead of 2.88 / 4.01). | `common.pss_asymptotic_F` |
| 4 | Long-run coefficients were taken from `ci_params` (theta / phi), which has the opposite sign to the long-run coefficient -theta / phi. | Signs in the long-run table were reversed. | `common.long_run` |
| 5 | Lag orders were compared by AIC on samples of different length. | Selected orders change for M1t, M3 and M5. | `econ.run_model` (common hold-back) |
| 6 | The Johansen system included the two I(0) variables. | Trace statistics rejected mechanically. On the I(1) variables there is no rejection, corrected or not. | `econ.py`, section 6 |
| 7 | The Gregory-Hansen statistic used an AIC-selected lag while its critical values were simulated with one lag, and the tabulated values could not be reproduced from the deposited script. The validation was against Gregory and Hansen (1996a), Table 1, not 1996b. | Statistic and critical values now come from one routine. No decision changes. | `gh_critical_values.py` |
| 8 | The text stated that unconstrained AIC selection prefers a pure autoregression. This is not reproduced. | The statement was removed from the manuscript. | - |

None of the corrections changes the paper's conclusion: no specification
supports a level relationship.

### Additions

* t-test of Pesaran, Shin and Smith (2001) on the lagged dependent level, with
  bounds simulated for the exact sample size (`common.simulate_bounds_cv`).
* Ten further specifications (seven listed here, M10 to M12 below): M4T (Tirana index), M6 (quarter-on-quarter
  credit growth), M7 (log levels, no common denominator), M8 (foreign-currency
  share), M9 (exposure shares only), M1r (regime dummy), M1s (seasonal dummies).
* M10, M11 and M12 add the standard macroeconomic determinants requested in review: the
  logarithm of the lek/euro exchange rate (`data/raw/boa_exchange_rates_end_of_month.csv`),
  the quarterly Labour Force Survey unemployment rate
  (`data/raw/instat_lfs_quarterly_indicators.xml`), and both together.
* Lag-order grid for M1, p = 1..4 and q = 1..3, on a common sample.
* Breusch-Godfrey test in chi-square and F form; OLS-based CUSUM test.
* `regimes.py`: interaction and Chow tests at 2022Q1 and 2020Q2, moratorium
  interaction, sup-F test with recursive bootstrap, sub-sample bounds tests.
* `size_power.py`: data-generating processes calibrated to the estimated
  models; designs with resampled residuals, AR(1) errors, a break in drift and
  simulated regressors; power over adjustment speed and long-run coefficient.
* `scenario.py`: new impairments separated from the provisioned stock;
  impairment rate conditional on negative equity; loan-to-value ratios
  distributed within buckets; workout cost and time to recovery; one-at-a-time
  sensitivity; Monte Carlo loss distribution; reverse stress test; benchmark
  under the IMF staff assumptions. The version 1 calculation is kept as
  `run_v1` for reconciliation.
* `build.py`: foreign-currency shares of real estate loans, logarithms of real
  estate and other loans, quarter-on-quarter credit growth, regime dummies, the
  lek/euro exchange rate and the unemployment rate.
  All version 1 columns are unchanged.

## 1.1.0 (first submission, October 2026)

Initial deposit.
