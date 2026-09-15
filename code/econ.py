import pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
RAW  = ROOT / "data" / "raw"
DER  = ROOT / "data" / "derived"
OUT  = ROOT / "output"
import warnings, json, math
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import statsmodels.api as sm
from statsmodels.tsa.stattools import adfuller
from statsmodels.tsa.ardl import ARDL, UECM, ardl_select_order
from statsmodels.stats.diagnostic import (acorr_breusch_godfrey, het_arch, linear_reset,
                                          breaks_cusumolsresid)
from statsmodels.stats.stattools import jarque_bera
from statsmodels.tsa.vector_ar.vecm import coint_johansen
from arch.unitroot import DFGLS, PhillipsPerron, KPSS, ZivotAndrews
from scipy import stats

R = {}
d = pd.read_csv(str(DER / "dataset_quarterly_2016Q1_2026Q2.csv"), index_col=0)
d.index = pd.PeriodIndex(d.index.str.replace("Q", "Q"), freq="Q")

# ------------------------------------------------------------------ 1. UNIT ROOTS
def ur_battery(s, name):
    s = s.dropna()
    out = {"n": len(s)}
    a = adfuller(s, autolag="AIC", regression="c");  out["adf"] = (a[0], a[1])
    g = DFGLS(s, trend="c");                          out["dfgls"] = (g.stat, g.pvalue)
    p = PhillipsPerron(s, trend="c");                 out["pp"] = (p.stat, p.pvalue)
    k = KPSS(s, trend="c");                           out["kpss"] = (k.stat, k.pvalue)
    try:
        z = ZivotAndrews(s, trend="c");               out["za"] = (z.stat, z.pvalue, str(s.index[z.break_idx] if hasattr(z, 'break_idx') else ""))
    except Exception:
        out["za"] = (np.nan, np.nan, "")
    return out

VARS = ["npl_logit", "ln_stok", "hip_share", "ppb_share", "eksp_ndertim", "ln_hpi",
        "pbb", "repo_real", "rritje_kredi", "car"]
ur = {}
for v in VARS:
    ur[v] = {"level": ur_battery(d[v], v), "diff": ur_battery(d[v].diff(), v)}
R["unit_roots"] = ur

# ------------------------------------------------------------------ 2. ARDL / BOUNDS
def run_ardl(dep, regs, sample, trend="c", maxlag=2, label=""):
    dd = d.loc[sample[0]:sample[1], [dep] + regs + ["d_covid"]].dropna()
    y = dd[dep]; X = dd[regs]; fx = dd[["d_covid"]]
    # (a) unrestricted AIC selection (may drop regressors) -> report
    sel = ardl_select_order(y, maxlag, X, maxlag, trend=trend, ic="aic", fixed=fx, glob=False)
    unres = (sel.model.ar_lags, {k: v for k, v in sel.model.dl_lags.items()})
    # (b) restricted grid: every regressor enters (lag 0..maxlag), ar 1..maxlag
    best = None
    for p in range(1, maxlag + 1):
        for q in range(1, maxlag + 1):
            mm = ARDL(y, p, X, q, trend=trend, fixed=fx)
            r = mm.fit()
            if best is None or r.aic < best[0]: best = (r.aic, p, q, mm, r)
    aic, p, q, mm, m = best
    u = UECM.from_ardl(mm).fit()
    case = 3 if trend == "c" else 5
    bt_a = u.bounds_test(case=case, asymptotic=True)
    bt_f = u.bounds_test(case=case, asymptotic=False, nsim=20000, rng=42)
    ols = sm.OLS(u.model.endog, u.model.exog).fit()
    res = ols.resid
    lm = acorr_breusch_godfrey(ols, nlags=4)
    arch = het_arch(res, nlags=4)
    jb = jarque_bera(res)
    try:
        rs = linear_reset(ols, power=2, use_f=True); reset = (float(rs.fvalue), float(rs.pvalue))
    except Exception:
        reset = (np.nan, np.nan)
    cus = breaks_cusumolsresid(ols.resid, ddof=ols.df_model)
    lr = {k: (float(u.ci_params[k]), float(u.ci_bse[k]), float(u.ci_pvalues[k])) for k in u.ci_params.index}
    ect_name = [c for c in u.params.index if c.startswith(dep + ".L1")]
    ect = (float(u.params[ect_name[0]]), float(u.bse[ect_name[0]]), float(u.pvalues[ect_name[0]])) if ect_name else None
    out = dict(label=label, dep=dep, regs=regs, sample=(str(dd.index[0]), str(dd.index[-1])), n=len(dd),
               unrestricted_aic_choice=str(unres), p=p, q=q, aic=float(aic), trend=trend,
               F=float(bt_a.statistic), p_asym=dict(bt_a.pvalue) if hasattr(bt_a.pvalue,'items') else str(bt_a.pvalue),
               p_fin=dict(bt_f.pvalue) if hasattr(bt_f.pvalue,'items') else str(bt_f.pvalue),
               cv_asym=bt_a.critical_values.to_dict() if hasattr(bt_a.critical_values, "to_dict") else str(bt_a.critical_values),
               cv_fin=bt_f.critical_values.to_dict() if hasattr(bt_f.critical_values, "to_dict") else str(bt_f.critical_values),
               ect=ect, longrun=lr,
               lm=(float(lm[0]), float(lm[1])), arch=(float(arch[0]), float(arch[1])),
               jb=(float(jb[0]), float(jb[1])), reset=reset, cusum=(float(cus[0]), float(cus[1])),
               r2=float(ols.rsquared), nparams=int(len(u.params)))
    return out, u, m

specs = [
 ("M1", "npl_logit", ["hip_share", "ppb_share", "pbb", "repo_real"], ("2016Q1", "2026Q1"), "c"),
 ("M1t", "npl_logit", ["hip_share", "ppb_share", "pbb", "repo_real"], ("2016Q1", "2026Q1"), "ct"),
 ("M2", "ln_stok", ["hip_share", "ppb_share", "pbb", "repo_real"], ("2016Q1", "2026Q1"), "c"),
 ("M3", "npl_logit", ["eksp_ndertim", "pbb", "repo_real"], ("2016Q3", "2026Q1"), "c"),
 ("M4", "npl_logit", ["ln_hpi", "pbb", "repo_real"], ("2020Q2", "2025Q4"), "c"),
 ("M5", "npl_logit", ["hip_share", "ppb_share", "pbb", "repo_real"], ("2017Q1", "2026Q1"), "c"),
]
R["ardl"] = {}
models = {}
for lab, dep, regs, smp, tr in specs:
    try:
        o, u, m = run_ardl(dep, regs, smp, trend=tr, label=lab); R["ardl"][lab] = o; models[lab] = (u, m)
        print(f"{lab}: n={o['n']} F={o['F']:.3f} ECT={o['ect']} p={o['p']} q={o['q']} | unrestricted AIC choice: {o['unrestricted_aic_choice']}")
        print("   p_asym:", o["p_asym"], "| p_fin:", o["p_fin"]); print("   fin-sample CV:", o["cv_fin"])
        print("   LR:", {k: (round(v[0],3), round(v[2],3)) for k, v in o["longrun"].items()})
        print("   diag LM p=%.3f ARCH p=%.3f JB p=%.3f RESET p=%.3f CUSUM p=%.3f R2=%.3f" % (o['lm'][1],o['arch'][1],o['jb'][1],o['reset'][1],o['cusum'][1],o['r2']))
    except Exception as e:
        print(lab, "FAILED", repr(e)); R["ardl"][lab] = {"error": repr(e)}

# ------------------------------------------------------------------ 3. GREGORY-HANSEN (level shift & regime shift)
def gh_test(y, X, trim=0.15):
    y = np.asarray(y); X = np.asarray(X); n = len(y)
    best = {"C": (np.inf, None), "CS": (np.inf, None)}
    for tb in range(int(trim*n), int((1-trim)*n)):
        D = (np.arange(n) >= tb).astype(float)
        # Model C: level shift
        Z = np.column_stack([np.ones(n), D, X])
        e = y - Z @ np.linalg.lstsq(Z, y, rcond=None)[0]
        t = adfuller(e, regression="n", autolag="AIC")[0]
        if t < best["C"][0]: best["C"] = (t, tb)
        # Model C/S: regime shift
        Z = np.column_stack([np.ones(n), D, X, X * D[:, None]])
        e = y - Z @ np.linalg.lstsq(Z, y, rcond=None)[0]
        t = adfuller(e, regression="n", autolag="AIC")[0]
        if t < best["CS"][0]: best["CS"] = (t, tb)
    return best
gh = {}
for lab, dep, regs, smp, tr in specs[:1] + specs[2:4]:
    dd = d.loc[smp[0]:smp[1], [dep] + regs].dropna()
    b = gh_test(dd[dep].values, dd[regs].values)
    gh[lab] = {k: (float(v[0]), str(dd.index[v[1]])) for k, v in b.items()}
    gh[lab]["m"] = len(regs)
    print("GH", lab, gh[lab])
R["gregory_hansen"] = gh
# Gregory-Hansen ADF* critical values, simulated for n=42 by Monte Carlo (05_gh_cv.py):
# (1%, 5%, 10%). Simulation validated against Gregory and Hansen (1996b, OBES Table 1).
R["gh_cv"] = {"C":  {3:(-6.125,-5.481,-5.143), 4:(-6.486,-5.748,-5.452)},
              "CS": {3:(-7.015,-6.160,-5.805), 4:(-7.488,-6.543,-6.190)}}
R["gh_cv_source"] = "Monte Carlo, 2000 replications, n=42, 15% trimming, one lag; validated against Gregory and Hansen (1996b) Table 1"

# ------------------------------------------------------------------ 4. JOHANSEN cross-check (M1 set)
dd = d.loc["2016Q1":"2026Q1", ["npl_logit","hip_share","ppb_share","pbb","repo_real"]].dropna()
jo = coint_johansen(dd.values, det_order=0, k_ar_diff=1)
R["johansen"] = {"trace": jo.lr1.tolist(), "cv_trace_95": jo.cvt[:,1].tolist(),
                 "maxeig": jo.lr2.tolist(), "cv_maxeig_95": jo.cvm[:,1].tolist(), "n": len(dd)}
print("Johansen trace:", np.round(jo.lr1,2), "cv95:", np.round(jo.cvt[:,1],2))

# ------------------------------------------------------------------ 5. PLACEBO: random-walk regressors
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
R["placebo"] = {"N": N, "n": n, "share_sig_levels": sig_lvl/N, "mean_r2_levels": float(np.mean(r2s)),
                "median_r2_levels": float(np.median(r2s)), "share_sig_diff": sig_dif/N}
print("Placebo:", R["placebo"])

# ------------------------------------------------------------------ 6. CORRELATIONS levels vs differences
V = ["hip_share","ppb_share","eksp_ndertim","ln_hpi","pbb","repo_real","rritje_kredi"]
R["corr"] = {v: (float(d["npl_logit"].corr(d[v])), float(d["npl_logit"].diff().corr(d[v].diff()))) for v in V}

# ------------------------------------------------------------------ 7. POWER / MDE (M1 long-run coefficients)
u = models["M1"][0]
dfree = int(u.df_resid)
tcrit = stats.t.ppf(0.975, dfree) + stats.t.ppf(0.80, dfree)
R["mde"] = {k: float(tcrit * v[1]) for k, v in R["ardl"]["M1"]["longrun"].items()}
R["mde_df"] = dfree
print("MDE:", R["mde"])

# ------------------------------------------------------------------ 8. Descriptives
desc = d[["npl","npl_logit","ln_stok","hip_share","ppb_share","eksp_ndertim","hpi_vendi","hpi_tirane","pbb","repo","inflacion","repo_real","car","npl_kap"]].describe().T[["count","mean","std","min","max"]]
R["desc"] = desc.round(3).to_dict(orient="index")

json.dump(R, open(str(OUT / "results.json"),"w"), indent=1, default=str)
print("saved")
