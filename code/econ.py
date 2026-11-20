"""
Estimation script, revision 2.

Changes relative to the version deposited with the original submission
(see docs/CHANGELOG.md for the full account):

  * Residual diagnostics are computed on the estimated error-correction
    model. Revision 1 computed them on `u.model.exog`, which is the static
    levels regression (erroneous Table 8).
  * The bounds F-statistic is computed on the model as specified, including
    the COVID dummy. Revision 1 used statsmodels' `bounds_test`, which drops
    fixed regressors.
  * Asymptotic critical values use k = number of regressors. Revision 1
    reported the values for k + 1.
  * The t-statistic on the lagged dependent level is compared with simulated
    PSS t-bounds, in addition to the F-test.
  * Lag orders are compared on a common estimation sample.
  * Johansen and Gregory-Hansen cross-checks use I(1) variables only.
  * Additional specifications requested in review (credit growth, a
    normalisation without the common denominator, euroisation, exposure-only,
    regime dummy, Tirana house prices) and a lag-order sensitivity grid.
"""
import json, math
import numpy as np, pandas as pd
import statsmodels.api as sm
from statsmodels.tsa.stattools import adfuller
from statsmodels.tsa.ardl import ARDL, UECM
from statsmodels.tsa.vector_ar.vecm import coint_johansen
from arch.unitroot import DFGLS, PhillipsPerron, KPSS, ZivotAndrews
from scipy import stats

from common import (load, fit_ecm, bounds_stats, long_run, diagnostics, bounds_cv,
                    pss_asymptotic_F, jdump, OUT)
from gh_critical_values import gh_stat

R = {}
d = load()

# ------------------------------------------------------------------ 1. UNIT ROOTS
def ur_battery(s):
    s = s.dropna()
    out = {"n": len(s)}
    a = adfuller(s, autolag="AIC", regression="c");  out["adf"] = (a[0], a[1])
    g = DFGLS(s, trend="c");                          out["dfgls"] = (g.stat, g.pvalue)
    p = PhillipsPerron(s, trend="c");                 out["pp"] = (p.stat, p.pvalue)
    k = KPSS(s, trend="c");                           out["kpss"] = (k.stat, k.pvalue)
    try:
        z = ZivotAndrews(s, trend="c");               out["za"] = (z.stat, z.pvalue)
    except Exception:
        out["za"] = (np.nan, np.nan)
    return out

VARS = ["npl_logit", "ln_stok", "hip_share", "ppb_share", "eksp_ndertim", "ln_hpi", "ln_hpi_tr",
        "pbb", "repo_real", "rritje_kredi", "dln_kredi", "re_fx_share", "ln_eur", "papunesia", "ln_re", "ln_other", "car"]
R["unit_roots"] = {v: {"level": ur_battery(d[v]), "diff": ur_battery(d[v].diff())} for v in VARS}

def integration_order(v):
    u = R["unit_roots"][v]
    if u["level"]["dfgls"][1] < 0.05 and u["level"]["kpss"][1] > 0.05:
        return "I(0)"
    return "I(1)"
R["integration"] = {v: integration_order(v) for v in VARS}
print("integration:", R["integration"])

# ------------------------------------------------------------------ 2. ARDL / BOUNDS
NSIM = 40000

def run_model(label, dep, regs, sample, trend="c", fixed=("d_covid",), pgrid=(1, 2), qgrid=(1, 2),
              note=""):
    cols = [dep] + list(regs) + list(fixed)
    dd = d.loc[sample[0]:sample[1], cols].dropna()
    hb = max(max(pgrid), max(qgrid))                 # common estimation sample for lag selection
    best = None
    for p in pgrid:
        for q in qgrid:
            res, y, X = fit_ecm(dd, dep, regs, p, q, trend=trend, fixed=fixed, hold_back=hb)
            if best is None or res.aic < best[0]:
                best = (res.aic, p, q, res, y, X)
    aic, p, q, res, y, X = best
    case = 3 if trend == "c" else 5
    k = len(regs)
    nobs = int(res.nobs)
    F, t, phi, se_phi, p_phi = bounds_stats(res, dep, regs)
    cv = bounds_cv(nobs, k, case, nsim=NSIM, seed=42)
    pF = (float((cv["_F_I0"] > F).mean()), float((cv["_F_I1"] > F).mean()))
    pt = (float((cv["_t_I0"] < t).mean()), float((cv["_t_I1"] < t).mean()))
    f_lo5, f_up5 = cv["F"]["5"]
    t_lo5, t_up5 = cv["t"]["5"]
    if F > f_up5 and t < t_up5:
        verdict = "Level relationship"
    elif F < f_lo5:
        verdict = "No"
    elif F > f_up5:
        verdict = "F rejects, t does not"
    else:
        verdict = "Inconclusive (F between bounds)"
    out = dict(label=label, dep=dep, regs=list(regs), fixed=list(fixed), note=note,
               sample=(str(dd.index[0]), str(dd.index[-1])), n=len(dd), nobs=nobs, k=k, case=case,
               p=p, q=q, aic=float(aic), bic=float(res.bic), trend=trend, nparams=len(res.params),
               F=F, t=t, phi=phi, se_phi=se_phi, p_phi=p_phi, pF=pF, pt=pt,
               cvF=cv["F"], cvt=cv["t"], cvF_asym=pss_asymptotic_F(k, case), verdict=verdict,
               longrun=long_run(res, dep, regs, trend=trend), diag=diagnostics(res),
               params={n: (float(b), float(s), float(pv)) for n, b, s, pv in
                       zip(res.names, res.params, res.bse, res.pvalues)})
    return out, res

M1REGS = ["hip_share", "ppb_share", "pbb", "repo_real"]
FULL = ("2016Q1", "2026Q1")
PRE = [
 ("M1",  "npl_logit", M1REGS, FULL, "c",  ("d_covid",), "H1, intercept (Case III)"),
 ("M1t", "npl_logit", M1REGS, FULL, "ct", ("d_covid",), "H1 with linear trend (Case V)"),
 ("M2",  "ln_stok",   M1REGS, FULL, "c",  ("d_covid",), "log NPL stock as dependent variable"),
 ("M3",  "npl_logit", ["eksp_ndertim", "pbb", "repo_real"], ("2016Q3", "2026Q1"), "c", ("d_covid",),
  "construction and real estate sector exposure"),
 ("M4",  "npl_logit", ["ln_hpi", "pbb", "repo_real"], ("2020Q2", "2025Q4"), "c", ("d_covid",),
  "national house price index (descriptive)"),
 ("M5",  "npl_logit", M1REGS, ("2017Q1", "2026Q1"), "c", ("d_covid",), "M1 on 2017Q1-2026Q1"),
]
ADD = [
 ("M4T", "npl_logit", ["ln_hpi_tr", "pbb", "repo_real"], ("2020Q2", "2025Q4"), "c", ("d_covid",),
  "Tirana house price index (descriptive)"),
 ("M6",  "npl_logit", M1REGS + ["dln_kredi"], FULL, "c", ("d_covid",),
  "M1 plus quarter-on-quarter credit growth"),
 ("M7",  "ln_stok", ["ln_re", "ln_other", "pbb", "repo_real"], FULL, "c", ("d_covid",),
  "log NPL stock on log real estate loans and log other loans (no common denominator)"),
 ("M8",  "npl_logit", M1REGS + ["re_fx_share"], FULL, "c", ("d_covid",),
  "M1 plus foreign-currency share of real estate loans"),
 ("M9",  "npl_logit", ["hip_share", "ppb_share"], FULL, "c", ("d_covid",),
  "exposure shares only (no year-on-year regressors)"),
 ("M10", "npl_logit", M1REGS + ["ln_eur"], FULL, "c", ("d_covid",),
  "M1 plus the log of the lek/euro exchange rate"),
 ("M11", "npl_logit", M1REGS + ["papunesia"], FULL, "c", ("d_covid",),
  "M1 plus the unemployment rate"),
 ("M12", "npl_logit", M1REGS + ["ln_eur", "papunesia"], FULL, "c", ("d_covid",),
  "M1 plus the exchange rate and the unemployment rate"),
 ("M1r", "npl_logit", M1REGS, FULL, "c", ("d_covid", "d_regime"),
  "M1 with a regime dummy from 2022Q1"),
 ("M1s", "npl_logit", M1REGS, FULL, "c", ("d_covid", "s1", "s2", "s3"),
  "M1 with quarterly seasonal dummies"),
]
R["ardl"] = {}
RES = {}
for lab, dep, regs, smp, tr, fx, note in PRE + ADD:
    o, res = run_model(lab, dep, regs, smp, trend=tr, fixed=fx, note=note)
    R["ardl"][lab] = o
    RES[lab] = res
    dg = o["diag"]
    print(f"{lab:4} n={o['n']} nobs={o['nobs']} ({o['p']},{o['q']}) F={o['F']:.3f} "
          f"[5%: {o['cvF']['5'][0]:.2f}/{o['cvF']['5'][1]:.2f}] t={o['t']:.2f} "
          f"[5%: {o['cvt']['5'][0]:.2f}/{o['cvt']['5'][1]:.2f}] phi={o['phi']:.3f} -> {o['verdict']}")
    print(f"      LM p={dg['lm'][1]:.3f} (F-form p={dg['lm_F'][1]:.3f}) ARCH p={dg['arch'][1]:.3f} JB p={dg['jb'][1]:.3f} "
          f"RESET p={dg['reset'][1]:.3f} CUSUM p={dg['cusum'][1]:.3f} adjR2={dg['r2_adj']:.3f}")

# check: the explicit design reproduces the statsmodels UECM estimates
o = R["ardl"]["M1"]
dd = d.loc[FULL[0]:FULL[1], ["npl_logit"] + M1REGS + ["d_covid"]].dropna()
u = UECM.from_ardl(ARDL(dd["npl_logit"], o["p"], dd[M1REGS], o["q"], trend="c", fixed=dd[["d_covid"]])).fit()
if int(u.nobs) == o["nobs"]:
    R["check_statsmodels_max_abs_diff"] = float(max(abs(o["params"][n][0] - u.params[n]) for n in u.params.index))
    print("max |coef difference| vs statsmodels UECM:", R["check_statsmodels_max_abs_diff"])

# what the submitted version reported, for the change log
ols_wrong = sm.OLS(u.model.endog, u.model.exog).fit()
R["v1_table8_explained"] = dict(
    regression_used="static OLS of the level of the dependent variable on the regressors in levels",
    nobs=int(ols_wrong.nobs), ncols=int(u.model.exog.shape[1]), ecm_nobs=o["nobs"], ecm_ncols=o["nparams"])

# ------------------------------------------------------------------ 3. LAG-ORDER SENSITIVITY (M1)
grid = []
dd = d.loc[FULL[0]:FULL[1], ["npl_logit"] + M1REGS + ["d_covid"]].dropna()
HB = 4
cvL = bounds_cv(len(dd) - HB, 4, 3, nsim=NSIM, seed=42)
for p in (1, 2, 3, 4):
    for q in (1, 2, 3):
        res, y, X = fit_ecm(dd, "npl_logit", M1REGS, p, q, trend="c", fixed=("d_covid",), hold_back=HB)
        F, t, phi, se, pv = bounds_stats(res, "npl_logit", M1REGS)
        grid.append(dict(p=p, q=q, nobs=int(res.nobs), nparams=len(res.params), dfr=int(res.df_resid),
                         aic=float(res.aic), bic=float(res.bic), F=F, t=t, phi=phi,
                         lm_p=diagnostics(res)["lm"][1], lmF_p=diagnostics(res)["lm_F"][1]))
R["lag_grid"] = dict(rows=grid, cvF=cvL["F"], cvt=cvL["t"], hold_back=HB)
ba = min(grid, key=lambda r: r["aic"]); bb = min(grid, key=lambda r: r["bic"])
R["lag_grid"]["aic_choice"] = (ba["p"], ba["q"]); R["lag_grid"]["bic_choice"] = (bb["p"], bb["q"])
R["lag_grid"]["max_F"] = max(r["F"] for r in grid)
print("lag grid: AIC", R["lag_grid"]["aic_choice"], "BIC", R["lag_grid"]["bic_choice"],
      "max F %.2f vs 5%% upper %.2f" % (R["lag_grid"]["max_F"], cvL["F"]["5"][1]))

# ------------------------------------------------------------------ 4. GREGORY-HANSEN, I(1) variables only
ghcv = json.load(open(str(OUT / "gh_cv.json")))["finite_sample"]
GH = [("M1", "npl_logit", ["hip_share", "ppb_share"], FULL),
      ("M2", "ln_stok", ["hip_share", "ppb_share"], FULL),
      ("M3", "npl_logit", ["eksp_ndertim"], ("2016Q3", "2026Q1"))]
gh = {}
for lab, dep, regs, smp in GH:
    dd = d.loc[smp[0]:smp[1], [dep] + regs].dropna()
    n, m = len(dd), len(regs)
    row = {"n": n, "m": m, "regs": regs}
    for model in ("C", "CS"):
        s, tb = gh_stat(dd[dep].values, dd[regs].values, model)
        row[model] = (s, str(dd.index[tb]))
        row[model + "_cv"] = ghcv[f"{model}_m{m}_n{n}"]
    gh[lab] = row
    print("GH", lab, {k: v for k, v in row.items() if k in ("C", "CS")}, row["C_cv"], row["CS_cv"])
R["gregory_hansen"] = gh

# ------------------------------------------------------------------ 5. JOHANSEN, I(1) variables only
def johansen(cols, smp, k_ar_diff=1):
    dd = d.loc[smp[0]:smp[1], cols].dropna()
    jo = coint_johansen(dd.values, det_order=0, k_ar_diff=k_ar_diff)
    T = len(dd) - (k_ar_diff + 1)                       # effective sample
    kk, pp = len(cols), k_ar_diff + 1
    ra = (T - kk * pp) / T
    return dict(vars=cols, n=len(dd), T_eff=T, lags_levels=pp, ra_factor=ra,
                trace=jo.lr1.tolist(), trace_ra=(jo.lr1 * ra).tolist(), cv_trace_95=jo.cvt[:, 1].tolist(),
                cv_trace_90=jo.cvt[:, 0].tolist(),
                maxeig=jo.lr2.tolist(), maxeig_ra=(jo.lr2 * ra).tolist(), cv_maxeig_95=jo.cvm[:, 1].tolist())
R["johansen"] = {
    "I1_only": johansen(["npl_logit", "hip_share", "ppb_share"], FULL, 1),
    "I1_only_2lags": johansen(["npl_logit", "hip_share", "ppb_share"], FULL, 2),
    "I1_only_stock": johansen(["ln_stok", "hip_share", "ppb_share"], FULL, 1),
    "v1_five_variables": johansen(["npl_logit"] + M1REGS, FULL, 1),
}
for k_, v in R["johansen"].items():
    print("Johansen", k_, "trace", np.round(v["trace"], 2), "RA", np.round(v["trace_ra"], 2),
          "cv95", np.round(v["cv_trace_95"], 2))

# ------------------------------------------------------------------ 6. PLACEBO: random-walk regressors
rng = np.random.default_rng(2026)
y = d["npl_logit"].dropna(); n = len(y); N = 5000
sig_lvl = 0; r2s = []; sig_dif = 0
dy = y.diff().dropna()
for i in range(N):
    rw = np.cumsum(rng.standard_normal(n))
    ols = sm.OLS(y.values, sm.add_constant(rw)).fit()
    sig_lvl += abs(ols.tvalues[1]) > 1.96; r2s.append(ols.rsquared)
    ols2 = sm.OLS(dy.values, sm.add_constant(np.diff(rw))).fit()
    sig_dif += abs(ols2.tvalues[1]) > 1.96
R["placebo"] = {"N": N, "n": n, "share_sig_levels": sig_lvl / N, "mean_r2_levels": float(np.mean(r2s)),
                "median_r2_levels": float(np.median(r2s)), "share_sig_diff": sig_dif / N}
print("Placebo:", R["placebo"])

# levels regression that corresponds to M1: the object revision 1 mistakenly diagnosed
dd = d.loc[FULL[0]:FULL[1], ["npl_logit"] + M1REGS].dropna()
lv = sm.OLS(dd["npl_logit"].values, sm.add_constant(dd[M1REGS].values)).fit()
from statsmodels.stats.stattools import durbin_watson
R["levels_regression_M1"] = dict(r2=float(lv.rsquared), dw=float(durbin_watson(lv.resid)), n=int(lv.nobs))

# ------------------------------------------------------------------ 6b. Exchange rate and unemployment: descriptive facts
R["eur_lek"] = {q: float(d.loc[q, "eur_lek"]) for q in ("2016Q1", "2019Q4", "2020Q1", "2021Q4", "2026Q1")}
R["papunesia"] = {q: float(d.loc[q, "papunesia"]) for q in ("2016Q1", "2021Q4", "2026Q1", "2026Q2")}
# The unemployment rate is not seasonally adjusted and rises in fourth quarters, when NPLs fall because of
# year-end write-offs. Correlation of quarterly changes before and after removing quarter-specific means:
_z = pd.DataFrame({"dn": d["npl_logit"].diff(), "du": d["papunesia"].diff()}).dropna()
_qtr = [p.quarter for p in _z.index]
_zs = _z - _z.groupby(_qtr).transform("mean")
R["un_seasonal"] = dict(corr_changes=float(_z["dn"].corr(_z["du"])), corr_changes_deseasonalised=float(_zs["dn"].corr(_zs["du"])),
                        mean_change_by_quarter={str(k): float(v) for k, v in _z["du"].groupby(_qtr).mean().items()})
print("unemployment, seasonal pattern:", R["un_seasonal"])

# ------------------------------------------------------------------ 7. CORRELATIONS levels vs differences
V = ["hip_share", "ppb_share", "eksp_ndertim", "ln_hpi", "ln_hpi_tr", "pbb", "repo_real", "rritje_kredi",
     "dln_kredi", "re_fx_share", "ln_eur", "papunesia"]
R["corr"] = {v: (float(d["npl_logit"].corr(d[v])), float(d["npl_logit"].diff().corr(d[v].diff()))) for v in V}

# ------------------------------------------------------------------ 8. HOUSE PRICES: semi-annual (non-interpolated) evidence
sem = d.loc[[i for i in d.index if i.quarter in (2, 4)]]
hp = {}
for col in ("ln_hpi", "ln_hpi_tr"):
    s = sem[[col, "npl_logit", "ln_stok"]].dropna()
    dz = s.diff().dropna()
    r = sm.OLS(dz["npl_logit"].values, sm.add_constant(dz[col].values)).fit()
    hp[col] = dict(n_levels=len(s), n_changes=len(dz),
                   corr_levels=float(s[col].corr(s["npl_logit"])),
                   corr_changes=float(dz[col].corr(dz["npl_logit"])),
                   corr_changes_stock=float(dz[col].corr(dz["ln_stok"])),
                   slope=float(r.params[1]), se=float(r.bse[1]), p=float(r.pvalues[1]),
                   n_falls=int((dz[col] < 0).sum()), min_change=float(dz[col].min()))
R["hpi_semiannual"] = hp
print("HPI semi-annual:", hp)

# ------------------------------------------------------------------ 9. MINIMUM DETECTABLE EFFECTS (M1)
o = R["ardl"]["M1"]
dfree = o["diag"]["dfr"]
tcrit = stats.t.ppf(0.975, dfree) + stats.t.ppf(0.80, dfree)
R["mde"] = {k_: float(tcrit * v[1]) for k_, v in o["longrun"].items()}
R["mde_df"] = dfree
print("MDE:", R["mde"], "df", dfree)

# ------------------------------------------------------------------ 10. Descriptives
desc = d[["npl", "npl_logit", "ln_stok", "hip_share", "ppb_share", "eksp_ndertim", "hpi_vendi", "hpi_tirane",
          "pbb", "repo", "inflacion", "repo_real", "car", "npl_kap", "dln_kredi", "hip_fx_share",
          "ppb_fx_share", "re_fx_share", "eur_lek", "papunesia"]].describe().T[["count", "mean", "std", "min", "max"]]
R["desc"] = desc.round(3).to_dict(orient="index")
R["fx_shares"] = {q: {c: float(d.loc[q, c]) for c in ("hip_fx_share", "ppb_fx_share", "re_fx_share")}
                  for q in ("2016Q1", "2021Q4", "2026Q2")}

dq = d["npl_logit"].diff(); ds = d["ln_stok"].diff()
R["seasonality"] = {"d_logit_by_quarter": {int(k_): float(v) for k_, v in dq.groupby(dq.index.quarter).mean().items()},
                    "d_lnstock_by_quarter": {int(k_): float(v) for k_, v in ds.groupby(ds.index.quarter).mean().items()}}

jdump(R, "results.json")
print("saved -> output/results.json")
