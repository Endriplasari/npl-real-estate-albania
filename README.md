# Real Estate Exposure and Nonperforming Loans in Albania

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22769974.svg)](https://doi.org/10.5281/zenodo.22769974)

Replication package for:

> Zeneli, M., Plasari, E., Meka, E., Baholli, I., & Kruja, S. *Real Estate Exposure
> and Nonperforming Loans in Albania: Limited Statistical Power and a Calibrated
> Stress Alternative.*

This repository contains everything needed to reproduce every number and table in
the paper, starting either from the raw statistical downloads or from the
assembled analysis dataset.

**Version 2.0.0** accompanies the revised manuscript. It corrects several
computations of version 1.1.0 and adds the analyses requested in review; see
`docs/CHANGELOG.md`. The outputs of the first submission are kept in
`output/v1_submitted/`.

---

## Quick start

```bash
git clone <repository-url>
cd <repository>
pip install -r requirements.txt
python code/run_all.py
```

Runtime is about four minutes on a standard machine, most of it the Gregory-Hansen
critical values. To reproduce only the
estimation results from the assembled dataset, skip `build.py` and run the other
scripts in the order given in `run_all.py`.

---

## Repository structure

```
code/
  build.py                 assembles the quarterly dataset from data/raw
  common.py                shared routines: conditional ECM, bounds F- and t-statistics,
                           simulated critical values, long-run coefficients, diagnostics
  gh_critical_values.py    Gregory-Hansen statistic and its Monte Carlo critical values
  econ.py                  unit roots, bounds tests, diagnostics, additional specifications,
                           lag grid, Gregory-Hansen, Johansen, placebo, minimum detectable effects
  regimes.py               structural stability across regimes and the COVID-19 period
  size_power.py            Monte Carlo size and power of the bounds procedure (calibrated)
  scenario.py              calibrated stress exercise, sensitivity, loss distribution,
                           reverse stress test, benchmark under IMF staff assumptions
  run_all.py               runs the above in order
data/
  raw/                     statistical downloads, unmodified
  derived/                 the assembled analysis dataset
output/                    generated results in JSON form
  v1_submitted/            outputs of version 1.1.0 (first submission)
docs/                      data provenance, house price reconstruction, change log
```

| Output file | Produced by | Tables in the revised paper |
|---|---|---|
| `output/results.json` | `econ.py` | 3, 7-14, 16, 19, A2; Section 6.12 |
| `output/gh_cv.json` | `gh_critical_values.py` | 13 (critical values) |
| `output/regimes.json` | `regimes.py` | 15 |
| `output/size_power.json` | `size_power.py` | 17, 18 |
| `output/scenario.json` | `scenario.py` | 4, 20-25 |

---

## Data

### What is included

`data/raw/` contains the statistical series exactly as downloaded, renamed for
legibility but otherwise unmodified. These are aggregate public statistics
published by the Bank of Albania and by INSTAT. They are reproduced here under
attribution so that the build step can be executed without manual retrieval.

`data/derived/dataset_quarterly_2016Q1_2026Q2.csv` is the analysis dataset:
42 quarterly observations, 2016Q1 to 2026Q2, 39 variables. It is produced by
`code/build.py` and is the only input to the estimation scripts.

### What is not included, and why

The Bank of Albania's *Results of the Survey on the Real Estate Market and the
House Price Index* is a copyrighted publication and is therefore linked rather
than redistributed. The survey reports publish percentage changes, not index
levels; the house price series used in the paper was reconstructed by chaining
those changes and anchoring them to the level reported against the 2013 base
period. The extracted figures, the issue from which each was taken and the
validation checks are documented in `docs/house_price_index.md`, and the
resulting series is embedded in `code/build.py`.

### Provenance

See `docs/data_provenance.md` for the source, the exact table name, the
retrieval path and the export settings for every series. One setting matters:
the Bank of Albania's time series interface defaults to zero decimal places,
which renders several series unusable. All exports here were taken with two
decimal places.

---

## Reproducibility

Random seeds are fixed throughout, so re-running the scripts reproduces the
published figures exactly:

| Seed       | Where                                        | What it controls                                   |
|------------|----------------------------------------------|----------------------------------------------------|
| `42`       | `common.simulate_bounds_cv`                  | finite-sample bounds F and t critical values       |
| `2026`     | `econ.py`                                    | random-walk placebo replications                   |
| `19960303` | `gh_critical_values.py`                      | Gregory-Hansen critical values                     |
| `20261114` | `regimes.py`, `size_power.py`, `scenario.py` | sup-F bootstrap, size and power, loss distribution |

Two sets of critical values are simulated rather than taken from published
tables. The bounds-test values follow the design of Pesaran, Shin and Smith
(2001) and the logic of Narayan (2005), and are generated for the exact number
of observations and regressors of each model (40,000 replications); at
T = 1,000 the routine reproduces the published asymptotic bounds. The
Gregory-Hansen values are generated for the sample sizes used (41 and 39
observations, 5,000 replications); the routine was validated by reproducing,
at n = 250, the asymptotic values in Gregory and Hansen (1996a, Table 1) for
m = 1 and m = 2 to within 0.22. The validation is stored in `output/gh_cv.json`.

Results were produced with Python 3.13, statsmodels 0.15.0 and arch 9.0.0.
Unit-root p-values can differ in the third decimal under other versions of
`arch`.

---

## Requirements

Python 3.11 or later. See `requirements.txt`.

---

## Citation

If you use this material, please cite the paper and this repository.
The archive is deposited on Zenodo under the concept DOI
[10.5281/zenodo.22769974](https://doi.org/10.5281/zenodo.22769974), which always
resolves to the most recent version. Author and version metadata are in
`CITATION.cff`.

## Licence

Code is released under the MIT Licence (`LICENSE`).
Data reproduced in `data/raw/` remain the property of their respective
publishers, the Bank of Albania and INSTAT, and are included under attribution
for the purpose of replication.
