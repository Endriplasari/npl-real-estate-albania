# Data provenance

All series are public. Retrieved September 2026.

## Bank of Albania — Statistics → Time Series database

Interface: Statistics → Time Series (`Seri_Kohore`). **Set decimal places to 2
before exporting.** The default of zero decimals renders the NPL ratio, the
capital adequacy ratio and all interest rate series unusable.

| File in `data/raw/` | Table in the database | Frequency | Coverage |
|---|---|---|---|
| `boa_financial_soundness_indicators.xlsx` | Financial soundness indicators (per cent) | Quarterly to Dec 2015, monthly thereafter | Mar 2014 – Jun 2026 |
| `boa_loans_by_purpose_currency.xlsx` | Loans by purpose of use and currency | Monthly | Dec 2015 – Jul 2026 |
| `boa_loans_by_district_activity.xlsx` | Loans by district and economic activity (NACE Rev.2) | Quarterly | 2016Q3 – 2026Q2 |
| `boa_loans_by_district_subject_maturity.xlsx` | Loans by district, subject and maturity | Quarterly | 2015Q4 – 2026Q2 |
| `boa_credit_to_economy.xlsx` | Credit to the economy | Monthly | Dec 2002 – Jul 2026 |
| `boa_cpi_2022_2026.xlsx` | Consumer price index (INSTAT series hosted by BoA) | Monthly | Jan 2022 – Aug 2026 |

Variables used from the financial soundness indicators table: item 2.6 (gross
non-performing loans to total loans), 1.1 (regulatory capital to risk-weighted
assets), 1.1.1 (regulatory capital), 1.1.2 (risk-weighted assets), 1.4.1 (NPLs
net of provisions to regulatory capital), 2.5 (return on assets).

Mortgage lending is the sum of items 1.3.1.4, 1.3.2.4, 1.3.3.4 and 1.3.4.4
(household loans for house purchase, by currency). Corporate real estate lending
is the sum of items 1.2.1.5, 1.2.2.5, 1.2.3.5 and 1.2.4.5.

Note: the financial soundness indicator series was revised in December 2017
following methodological alignment with the IMF compilation guide. The revision
concerns income-based indicators; series 2.6 appears unaffected but the break is
recorded in the paper.

## Bank of Albania — other

| File | Source | Note |
|---|---|---|
| `boa_policy_rates.xlsx` | Official rates of the Bank of Albania | A table of Supervisory Council decisions, not a time series. Converted to a monthly series by forward-fill in `build.py`. |

Interest rate statistics on new lending are **not** used. A regulation adopted in
September 2017 changed the compilation methodology, with first reporting in
December 2017 and a consolidation period to November 2018, producing an
unavoidable break. The policy rate is continuous over the sample and is used
instead.

## INSTAT

| File | Source | Frequency | Coverage |
|---|---|---|---|
| `instat_gdp_quarterly_production.xlsx` | Quarterly national accounts, production approach, real growth rates | Quarterly | 2000Q1 – 2026Q1 |
| `instat_cpi_annual_change.xlsx` | Consumer price index, annual change by COICOP group | Monthly | Jan 2017 – Dec 2025 |
| `instat_construction_cost_index.xlsx` | Construction cost index for new residential buildings | Quarterly | 2018Q1 – 2026Q2 |

Inflation for calendar 2016 is not covered by the file above. It is taken from
the INSTAT monthly releases *Indeksi i Çmimeve të Konsumit*, specifically the
November 2016 release (Table 2, annual changes for January to November) and the
December 2016 release (2.2 per cent). These values are hard-coded in `build.py`
with a comment recording the source.

Note: the CPI has been compiled with a new basket since January 2016
(December 2015 = 100, 333 items, ECOICOP classification).

## Not redistributed

| Source | Reason | Where to obtain |
|---|---|---|
| Survey on the Real Estate Market and the House Price Index, semi-annual issues 2020H1 to 2025H2 | Copyrighted publication | Bank of Albania, Financial Stability section and Monetary Policy → Surveys |
| Annual Supervision Report | Copyrighted publication | Bank of Albania |

The figures extracted from the survey issues are listed in
`docs/house_price_index.md` together with the issue each was taken from.
