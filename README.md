# Real Estate Exposure and Non-Performing Loans in Albania

Replication package for:

> Meka, E., Plasari, E., Zeneli, M., Baholli, I., & Kruja, S. *Real Estate Exposure
> and Non-Performing Loans in Albania: Identification Failure and a Calibrated
> Stress Alternative.*

This repository contains everything needed to reproduce every number, table and
figure in the paper, starting either from the raw statistical downloads or from
the assembled analysis dataset.

---

## Quick start

```bash
git clone <repository-url>
cd <repository>
pip install -r requirements.txt
python code/run_all.py
```

Runtime is roughly 25 minutes on a standard machine, dominated by the Monte Carlo
simulations. To reproduce only the estimation results from the assembled dataset,
skip `build.py` and run `econ.py` and `scenario.py` directly.

---

## Repository structure

```
code/
  build.py                 assembles the quarterly dataset from data/raw
  econ.py                  unit roots, bounds tests, cointegration, placebo, power
  gh_critical_values.py    Monte Carlo critical values for the Gregory-Hansen test
  scenario.py              calibrated collateral-coverage stress exercise
  run_all.py               runs the above in order
data/
  raw/                     statistical downloads, unmodified
  derived/                 the assembled analysis dataset
output/                    generated results in JSON form
docs/                      data provenance and the house price reconstruction
```

---

## Data

### What is included

`data/raw/` contains the statistical series exactly as downloaded, renamed for
legibility but otherwise unmodified. These are aggregate public statistics
published by the Bank of Albania and by INSTAT. They are reproduced here under
attribution so that the build step can be executed without manual retrieval.

`data/derived/dataset_quarterly_2016Q1_2026Q2.csv` is the analysis dataset:
42 quarterly observations, 2016Q1 to 2026Q2, 27 variables. It is produced by
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

Random seeds are fixed throughout: 20260909 for the bounds-test critical values,
19960303 for the Gregory-Hansen critical values, 4242 for the power analysis and
1000 onwards for the placebo replications. Re-running the scripts reproduces the
published figures exactly.

Two sets of critical values are simulated rather than taken from published
tables. The bounds-test values follow the logic of Narayan (2005) but are
generated for the exact sample length and regressor count used here. The
Gregory-Hansen values are generated for n = 42; the simulation code was
validated by reproducing, at n = 250, the asymptotic values published in
Gregory and Hansen (1996, Table 1) to within 0.18 across all cells tested.
That validation is in `output/gh_validation_CS.json`.

---

## Requirements

Python 3.11 or later. See `requirements.txt`.

---

## Citation

If you use this material, please cite the paper and this repository.
A DOI is minted for each release through Zenodo; see `CITATION.cff`.

## Licence

Code is released under the MIT Licence (`LICENSE`).
Data reproduced in `data/raw/` remain the property of their respective
publishers, the Bank of Albania and INSTAT, and are included under attribution
for the purpose of replication.
