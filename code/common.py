"""
Shared routines for the revised estimation code (revision 2).

Everything that the paper reports about the conditional error-correction
model (ECM) is computed here from an explicitly constructed design matrix,
so that the bounds F-statistic, the t-statistic on the lagged level of the
dependent variable, the long-run coefficients and the residual diagnostics
all refer to one and the same regression.

Why this file exists
--------------------
Revision 1 relied on statsmodels' UECM object in three places where its
behaviour differs from what the paper states:

  1. `u.model.exog` holds the original regressors in levels, not the ECM
     design. Diagnostics computed from it describe a static levels
     regression. (This produced the erroneous Table 8 of the submitted
     manuscript.)
  2. `UECMResults.bounds_test` re-estimates the model without the fixed
     regressors, so the reported F-statistic excluded the COVID dummy.
  3. `bounds_test(asymptotic=True)` looks up the tabulated critical values
     with k equal to the number of regressors plus one.

The functions below avoid all three. The design reproduces the statsmodels
UECM parameter estimates exactly (checked in econ.py), which shows that the
model is the same and only the auxiliary calculations differ.
"""
import pathlib, json, warnings
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.diagnostic import (acorr_breusch_godfrey, het_arch,
                                          linear_reset, breaks_cusumolsresid)
from statsmodels.stats.stattools import jarque_bera

warnings.filterwarnings("ignore")

ROOT = pathlib.Path(__file__).resolve().parents[1]
DER = ROOT / "data" / "derived"
OUT = ROOT / "output"

REGIME_START = "2022Q1"       # first quarter of the credit-expansion regime


def load():
    d = pd.read_csv(str(DER / "dataset_quarterly_2016Q1_2026Q2.csv"), index_col=0)
    d.index = pd.PeriodIndex(d.index, freq="Q")
    for k in (1, 2, 3):                                   # quarterly seasonal dummies (Q4 is the base)
        d[f"s{k}"] = (d.index.quarter == k).astype(float)
    return d


# --------------------------------------------------------------------------
# Conditional ECM
# --------------------------------------------------------------------------
def ecm_design(dd, dep, regs, p, q, trend="c", fixed=(), hold_back=None):
    """Design of the unrestricted conditional ECM implied by an ARDL(p, q,...,q).

        D.y_t = c (+ d*t) + phi*y_{t-1} + sum_j theta_j*x_{j,t-1}
                + sum_{i=1}^{p-1} g_i*D.y_{t-i}
                + sum_j sum_{i=0}^{q-1} h_{ji}*D.x_{j,t-i}
                + fixed_t + e_t

    `hold_back` fixes the number of initial observations that are dropped, so
    that models with different lag orders can be compared on one sample.
    Returns (y, X) as a Series and a DataFrame with named columns.
    """
    hb = max(p, q) if hold_back is None else hold_back
    assert hb >= max(p, q)
    y = dd[dep]
    dy = y.diff()
    cols = {}
    cols["const"] = pd.Series(1.0, index=dd.index)
    if trend == "ct":
        cols["trend"] = pd.Series(np.arange(1, len(dd) + 1, dtype=float), index=dd.index)
    cols[f"{dep}.L1"] = y.shift(1)
    for r in regs:
        cols[f"{r}.L1"] = dd[r].shift(1)
    for i in range(1, p):
        cols[f"D.{dep}.L{i}"] = dy.shift(i)
    for r in regs:
        dx = dd[r].diff()
        for i in range(0, q):
            cols[f"D.{r}.L{i}"] = dx.shift(i)
    for f in fixed:
        cols[f] = dd[f]
    X = pd.DataFrame(cols)
    yy = dy.rename(f"D.{dep}")
    X = X.iloc[hb:]
    yy = yy.iloc[hb:]
    return yy, X


def fit_ecm(dd, dep, regs, p, q, trend="c", fixed=(), hold_back=None):
    y, X = ecm_design(dd, dep, regs, p, q, trend=trend, fixed=fixed, hold_back=hold_back)
    res = sm.OLS(y.values, X.values).fit()
    res.names = list(X.columns)
    return res, y, X


def level_terms(names, dep, regs):
    return [names.index(f"{dep}.L1")] + [names.index(f"{r}.L1") for r in regs]


def bounds_stats(res, dep, regs):
    """Bounds F (all lagged levels jointly zero) and t (lagged dependent level)."""
    names = res.names
    lev = level_terms(names, dep, regs)
    R = np.zeros((len(lev), len(names)))
    for i, j in enumerate(lev):
        R[i, j] = 1.0
    F = float(np.squeeze(res.f_test(R).fvalue))
    j = names.index(f"{dep}.L1")
    return F, float(res.tvalues[j]), float(res.params[j]), float(res.bse[j]), float(res.pvalues[j])


def long_run(res, dep, regs, trend="c"):
    """Long-run coefficients -theta_j/phi with delta-method standard errors."""
    names = res.names
    jphi = names.index(f"{dep}.L1")
    phi = res.params[jphi]
    V = res.cov_params()
    out = {}
    terms = list(regs) + (["trend"] if trend == "ct" else [])
    for r in terms:
        jr = names.index(f"{r}.L1") if r != "trend" else names.index("trend")
        th = res.params[jr]
        b = -th / phi
        g = np.zeros(len(names))
        g[jr] = -1.0 / phi
        g[jphi] = th / phi ** 2
        se = float(np.sqrt(g @ V @ g))
        z = b / se
        from scipy import stats as st
        out[r] = (float(b), se, float(2 * (1 - st.norm.cdf(abs(z)))))
    return out


def diagnostics(res):
    """Residual diagnostics of the estimated ECM itself.

    The Breusch-Godfrey test is reported in its chi-square (LM) form and in its
    F form; the latter is better behaved when the number of parameters is large
    relative to the sample (Kiviet, 1986).
    """
    e = res.resid
    lm = acorr_breusch_godfrey(res, nlags=4)
    arch = het_arch(e, nlags=4)
    jb = jarque_bera(e)
    try:
        rs = linear_reset(res, power=2, use_f=True)
        reset = (float(rs.fvalue), float(rs.pvalue))
    except Exception:
        reset = (float("nan"), float("nan"))
    cus = breaks_cusumolsresid(e, ddof=int(res.df_model) + 1)
    return dict(lm=(float(lm[0]), float(lm[1])), lm_F=(float(lm[2]), float(lm[3])), arch=(float(arch[0]), float(arch[1])),
                jb=(float(jb[0]), float(jb[1])), reset=reset,
                cusum=(float(cus[0]), float(cus[1])),
                r2=float(res.rsquared), r2_adj=float(res.rsquared_adj),
                sigma=float(np.sqrt(res.scale)), dfr=int(res.df_resid))


# --------------------------------------------------------------------------
# Critical values for the bounds F- and t-tests, by simulation
# --------------------------------------------------------------------------
def simulate_bounds_cv(nobs, k, case, nsim=40000, seed=42):
    """Finite-sample critical values of the PSS bounds F- and t-statistics.

    The data-generating process and the test regression follow Pesaran, Shin
    and Smith (2001): y is a driftless random walk, the k regressors are
    independent random walks (upper bound, I(1)) or white noise (lower bound,
    I(0)), and D.y_t is regressed on the deterministic terms, y_{t-1} and
    x_{t-1}. `nobs` is the number of observations entering the regression.
    case 3: unrestricted intercept; case 5: unrestricted intercept and trend.
    """
    rs = np.random.default_rng(seed)
    ndet = 1 if case == 3 else 2
    block = max(1, 60_000_000 // (8 * nobs * (k + 1 + ndet)))
    F = {"I0": np.empty(nsim), "I1": np.empty(nsim)}
    T = {"I0": np.empty(nsim), "I1": np.empty(nsim)}
    done = 0
    while done < nsim:
        b = min(block, nsim - done)
        e = rs.standard_normal((b, nobs + 1, k + 1))
        y = np.cumsum(e[:, :, 0], axis=1)
        lhs = np.diff(y, axis=1)[:, :, None]
        for kind in ("I1", "I0"):
            x = np.cumsum(e[:, :, 1:], axis=1) if kind == "I1" else e[:, :, 1:]
            rhs = np.empty((b, nobs, k + 1 + ndet))
            rhs[:, :, 0] = y[:, :-1]
            rhs[:, :, 1:k + 1] = x[:, :-1]
            rhs[:, :, k + 1] = 1.0
            if ndet == 2:
                rhs[:, :, k + 2] = np.arange(nobs, dtype=float)
            rt = np.transpose(rhs, (0, 2, 1))
            xpx = rt @ rhs
            xpxi = np.linalg.inv(xpx)
            beta = xpxi @ (rt @ lhs)
            u = lhs - rhs @ beta
            rss_u = (u ** 2).sum(axis=(1, 2))
            dfree = nobs - (k + 1 + ndet)
            s2 = rss_u / dfree
            det = rhs[:, :, k + 1:]
            dt = np.transpose(det, (0, 2, 1))
            br = np.linalg.solve(dt @ det, dt @ lhs)
            ur = lhs - det @ br
            rss_r = (ur ** 2).sum(axis=(1, 2))
            F[kind][done:done + b] = ((rss_r - rss_u) / (k + 1)) / s2
            T[kind][done:done + b] = beta[:, 0, 0] / np.sqrt(s2 * xpxi[:, 0, 0])
        done += b
    out = {"nobs": nobs, "k": k, "case": case, "nsim": nsim, "seed": seed, "F": {}, "t": {}}
    for lvl in (10, 5, 1):
        out["F"][str(lvl)] = (float(np.percentile(F["I0"], 100 - lvl)), float(np.percentile(F["I1"], 100 - lvl)))
        out["t"][str(lvl)] = (float(np.percentile(T["I0"], lvl)), float(np.percentile(T["I1"], lvl)))
    out["_F_I1"] = F["I1"]
    out["_F_I0"] = F["I0"]
    out["_t_I1"] = T["I1"]
    out["_t_I0"] = T["I0"]
    return out


_CV_CACHE = {}


def bounds_cv(nobs, k, case, nsim=40000, seed=42):
    key = (nobs, k, case, nsim, seed)
    if key not in _CV_CACHE:
        _CV_CACHE[key] = simulate_bounds_cv(nobs, k, case, nsim=nsim, seed=seed)
    return _CV_CACHE[key]


def pss_asymptotic_F(k, case):
    """Asymptotic PSS F-bounds, indexed by the number of regressors k."""
    from statsmodels.tsa.ardl import pss_critical_values as pc
    pct = list(pc.crit_percentiles)
    lo = pc.crit_vals[(k, case, False)]
    up = pc.crit_vals[(k, case, True)]
    return {str(100 - p): (float(lo[i]), float(up[i])) for i, p in enumerate(pct) if p in (90, 95, 99)}


def jdump(obj, name):
    def clean(o):
        if isinstance(o, dict):
            return {str(k): clean(v) for k, v in o.items() if not str(k).startswith("_")}
        if isinstance(o, (list, tuple)):
            return [clean(v) for v in o]
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, (np.bool_,)):
            return bool(o)
        return o
    json.dump(clean(obj), open(str(OUT / name), "w"), indent=1, default=str)
