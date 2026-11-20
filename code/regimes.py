"""
Structural stability of the exposure-NPL relation (revision 2).

Added in response to the review:
  * the 2021/2022 change from a write-off regime to a credit-expansion regime
    is modelled explicitly (interaction terms, Chow tests, sub-sample bounds
    tests);
  * the COVID-19 period is treated as a candidate break (interaction of real
    estate exposure with the pandemic period and with the moratorium quarters);
  * a sup-F test with an unknown break date, with a bootstrap p-value.

The analysis is run on a first-difference specification, whose validity does
not depend on the existence of a level relationship:

    D.y_t = a + b*D.exposure_t + c1*gdp_t + c2*realrate_t + g*covid_t
            + r*D.y_{t-1} + e_t

where exposure is the total real estate share of the loan book (mortgage plus
business real estate share).
"""
import numpy as np, pandas as pd
import statsmodels.api as sm
from scipy import stats

from common import load, fit_ecm, bounds_stats, bounds_cv, diagnostics, jdump, REGIME_START

SEED = 20261114
NBOOT = 2000
S = {}
d = load()

x = d.loc["2016Q1":"2026Q1"].copy()
for dep in ("npl_logit", "ln_stok"):
    x[f"d_{dep}"] = x[dep].diff()
    x[f"d_{dep}_1"] = x[f"d_{dep}"].shift(1)
x["dexp"] = x["eksp_total"].diff()
x["const"] = 1.0
BASE = ["const", "dexp", "pbb", "repo_real", "d_covid"]


def fd_model(dep, D=None, interact=("const", "dexp")):
    cols = BASE + [f"d_{dep}_1"]
    z = x[[f"d_{dep}"] + cols + (["d_regime", "d_postcovid"])].dropna()
    X = z[cols].copy()
    extra = []
    if D is not None:
        Dv = z[D] if isinstance(D, str) else D.reindex(z.index)
        for c in interact:
            name = f"D*{c}"
            X[name] = Dv.values * z[c].values
            extra.append(name)
    y = z[f"d_{dep}"]
    ols = sm.OLS(y.values, X.values).fit()
    hac = sm.OLS(y.values, X.values).fit(cov_type="HAC", cov_kwds={"maxlags": 3})
    names = list(X.columns)
    out = dict(n=int(ols.nobs), names=names,
               coef={n_: (float(b), float(s), float(p), float(sh), float(ph))
                     for n_, b, s, p, sh, ph in zip(names, ols.params, ols.bse, ols.pvalues, hac.bse, hac.pvalues)},
               r2=float(ols.rsquared), sigma=float(np.sqrt(ols.scale)))
    if extra:
        Rm = np.zeros((len(extra), len(names)))
        for i, e in enumerate(extra):
            Rm[i, names.index(e)] = 1.0
        ft = ols.f_test(Rm)
        out["chow_F"] = float(np.squeeze(ft.fvalue)); out["chow_p"] = float(ft.pvalue)
        out["chow_df"] = (len(extra), int(ols.df_resid))
        w = hac.wald_test(Rm, use_f=True, scalar=True)
        out["chow_F_hac"] = float(w.statistic); out["chow_p_hac"] = float(w.pvalue)
    return out, ols, z, X


# ---------------------------------------------------------------- A. interaction / Chow tests
A = {}
for dep in ("npl_logit", "ln_stok"):
    base, ols0, z0, X0 = fd_model(dep)
    A[dep] = {"baseline": base}
    for lab, D in (("regime_2022Q1", "d_regime"), ("covid_onset_2020Q2", "d_postcovid")):
        o, _, _, _ = fd_model(dep, D=D)
        A[dep][lab] = o
    # moratorium quarters only: slope interaction (the level effect is already in the baseline)
    o, _, _, _ = fd_model(dep, D="d_covid", interact=("dexp",))
    A[dep]["moratorium_quarters"] = o
    for lab in ("regime_2022Q1", "covid_onset_2020Q2", "moratorium_quarters"):
        r = A[dep][lab]
        print(f"{dep:10} {lab:22} Chow F={r['chow_F']:.2f} p={r['chow_p']:.3f} (HAC p={r['chow_p_hac']:.3f}) | "
              + "  ".join(f"{k}={v[0]:.3f} (p={v[2]:.3f})" for k, v in r["coef"].items() if k.startswith("D*")))
S["first_difference"] = A

# ---------------------------------------------------------------- B. sup-F with unknown break date
def supF(yv, Xb, dates, cand):
    best, when = -np.inf, None
    k0 = Xb.shape[1]
    rss0 = sm.OLS(yv, Xb).fit().ssr
    path = []
    for tb in cand:
        D = (np.arange(len(yv)) >= tb).astype(float)
        Xa = np.column_stack([Xb, D, D * Xb[:, 1]])
        r = sm.OLS(yv, Xa).fit()
        F = ((rss0 - r.ssr) / 2) / (r.ssr / (len(yv) - Xa.shape[1]))
        path.append(float(F))
        if F > best:
            best, when = F, tb
    return float(best), when, path

SUP = {}
rng = np.random.default_rng(SEED)
for dep in ("npl_logit", "ln_stok"):
    cols = BASE + [f"d_{dep}_1"]
    z = x[[f"d_{dep}"] + cols].dropna()
    yv = z[f"d_{dep}"].values
    Xb = z[cols].values
    n = len(yv)
    cand = list(range(int(0.15 * n), int(0.85 * n)))
    F, when, path = supF(yv, Xb, z.index, cand)
    # recursive residual bootstrap under the null of no break
    fit0 = sm.OLS(yv, Xb).fit()
    b0, e0 = fit0.params, fit0.resid
    jlag = cols.index(f"d_{dep}_1")
    y_init = float(Xb[0, jlag])                        # D.y_{t-1} at the first observation
    boot = np.empty(NBOOT)
    for bi in range(NBOOT):
        es = rng.choice(e0, size=n, replace=True)
        ys = np.empty(n); Xs = Xb.copy(); prev = y_init
        for t in range(n):
            Xs[t, jlag] = prev
            ys[t] = Xs[t] @ b0 + es[t]
            prev = ys[t]
        boot[bi] = supF(ys, Xs, None, cand)[0]
    SUP[dep] = dict(supF=F, break_date=str(z.index[when]), n=n, trimming=0.15,
                    first_candidate=str(z.index[cand[0]]), last_candidate=str(z.index[cand[-1]]),
                    boot_p=float((boot >= F).mean()), boot_cv95=float(np.quantile(boot, 0.95)),
                    boot_cv90=float(np.quantile(boot, 0.90)), nboot=NBOOT,
                    path={str(z.index[c]): f for c, f in zip(cand, path)})
    print(f"sup-F {dep}: {F:.2f} at {SUP[dep]['break_date']}, bootstrap p={SUP[dep]['boot_p']:.3f}, "
          f"cv95={SUP[dep]['boot_cv95']:.2f}")
S["supF"] = SUP

# ---------------------------------------------------------------- C. sub-sample bounds tests (parsimonious)
SUB = {}
specs = [("write-off regime", "2016Q1", "2021Q4", ("d_covid",)),
         ("credit-expansion regime", REGIME_START, "2026Q1", ())]
for dep in ("npl_logit", "ln_stok"):
    for lab, a, b, fixed in specs:
        dd = d.loc[a:b, [dep, "eksp_total", "pbb"] + list(fixed)].dropna()
        res, y, X = fit_ecm(dd, dep, ["eksp_total", "pbb"], 1, 1, trend="c", fixed=fixed)
        F, t, phi, se, pv = bounds_stats(res, dep, ["eksp_total", "pbb"])
        cv = bounds_cv(int(res.nobs), 2, 3, nsim=40000, seed=42)
        dg = diagnostics(res)
        SUB[f"{dep}|{lab}"] = dict(dep=dep, regime=lab, sample=(a, b), n=len(dd), nobs=int(res.nobs),
                                   nparams=len(res.params), F=F, t=t, phi=phi, se_phi=se, cvF=cv["F"],
                                   cvt=cv["t"], lm_p=dg["lm"][1], cusum_p=dg["cusum"][1])
        print(f"sub-sample {dep:10} {lab:24} nobs={int(res.nobs)} F={F:.2f} [5%: {cv['F']['5'][0]:.2f}/"
              f"{cv['F']['5'][1]:.2f}] t={t:.2f} [5%: {cv['t']['5'][0]:.2f}/{cv['t']['5'][1]:.2f}] phi={phi:.3f}")
S["subsample_bounds"] = SUB

# ---------------------------------------------------------------- D. descriptive contrast of the two regimes
DES = {}
for lab, a, b in (("write-off regime", "2016Q2", "2021Q4"), ("credit-expansion regime", REGIME_START, "2026Q1")):
    z = x.loc[a:b]
    DES[lab] = dict(n=int(z["d_npl_logit"].notna().sum()),
                    mean_d_logit=float(z["d_npl_logit"].mean()), sd_d_logit=float(z["d_npl_logit"].std()),
                    mean_d_lnstock=float(z["d_ln_stok"].mean()), sd_d_lnstock=float(z["d_ln_stok"].std()),
                    mean_dln_credit=float(z["dln_kredi"].mean()),
                    mean_dexp=float(z["dexp"].mean()),
                    corr_dlogit_dexp=float(z["d_npl_logit"].corr(z["dexp"])),
                    corr_dlnstock_dexp=float(z["d_ln_stok"].corr(z["dexp"])))
    print(lab, {k: round(v, 3) for k, v in DES[lab].items()})
S["regime_descriptives"] = DES

jdump(S, "regimes.json")
print("saved -> output/regimes.json")
