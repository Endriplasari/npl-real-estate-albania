"""
Size and power of the ARDL bounds procedure under the M1 design (revision 2).

Revision 1 used idealised data-generating processes (independent unit-variance
random walks, independent errors). In response to the review, the experiments
are now calibrated to the estimated processes:

  * the null model for the dependent variable is the autoregression in first
    differences estimated on the data, with its estimated innovation variance;
  * the regressors are the observed ones (held fixed) or are simulated from
    processes estimated on them;
  * serially correlated errors and a break in the drift of the dependent
    variable at the 2022Q1 regime change are added as separate designs;
  * under the alternative, the equilibrium error has the residual standard
    deviation of the estimated error-correction model, and the long-run
    coefficient is stated in economically interpretable units.

The test applied in every replication is the one applied in the paper: the
conditional ECM with the lag orders selected for M1, an unrestricted intercept
and the COVID dummy, evaluated against the finite-sample 5% bounds for the
F-statistic and for the t-statistic on the lagged dependent level.
"""
import json
import numpy as np
import statsmodels.api as sm

from common import load, fit_ecm, bounds_stats, bounds_cv, diagnostics, jdump, OUT, REGIME_START

SEED = 20261114
NREP_SIZE = 5000
NREP_POWER = 2000
DEP = "npl_logit"
REGS = ["hip_share", "ppb_share", "pbb", "repo_real"]
I1 = [0, 1]            # columns of REGS treated as I(1)
I0 = [2, 3]            # columns of REGS treated as I(0)

rng = np.random.default_rng(SEED)
d = load()
dd = d.loc["2016Q1":"2026Q1", [DEP] + REGS + ["d_covid", "d_regime"]].dropna()
m1 = json.load(open(str(OUT / "results.json")))["ardl"]["M1"]
P, Q = m1["p"], m1["q"]
HB = 2
Y = dd[DEP].values.astype(float)
X = dd[REGS].values.astype(float)
COVID = dd["d_covid"].values.astype(float)
REGIME = dd["d_regime"].values.astype(float)
T = len(Y)


def ecm_test(y, Xm, dummy, p=P, q=Q, hb=HB):
    """F and t of the bounds procedure in the conditional ECM (intercept, dummy)."""
    n = len(y)
    dy = np.diff(y, prepend=np.nan)
    dX = np.diff(Xm, axis=0, prepend=np.full((1, Xm.shape[1]), np.nan))
    idx = np.arange(hb, n)
    cols = [np.ones(len(idx)), y[idx - 1]]
    for j in range(Xm.shape[1]):
        cols.append(Xm[idx - 1, j])
    nlev = 1 + Xm.shape[1]
    for i in range(1, p):
        cols.append(dy[idx - i])
    for j in range(Xm.shape[1]):
        for i in range(0, q):
            cols.append(dX[idx - i, j])
    cols.append(dummy[idx])
    Z = np.column_stack(cols)
    yv = dy[idx]
    ZtZi = np.linalg.inv(Z.T @ Z)
    b = ZtZi @ (Z.T @ yv)
    e = yv - Z @ b
    dfr = len(yv) - Z.shape[1]
    s2 = e @ e / dfr
    Rm = np.zeros((nlev, Z.shape[1]))
    for i in range(nlev):
        Rm[i, 1 + i] = 1.0
    rb = Rm @ b
    F = float(rb @ np.linalg.solve(Rm @ ZtZi @ Rm.T, rb) / nlev / s2)
    t = float(b[1] / np.sqrt(s2 * ZtZi[1, 1]))
    return F, t, dfr


# ------------------------------------------------------------ check against the estimation code
res, y_, X_ = fit_ecm(dd, DEP, REGS, P, Q, trend="c", fixed=("d_covid",), hold_back=HB)
F_ref, t_ref, phi_ref, _, _ = bounds_stats(res, DEP, REGS)
F_chk, t_chk, dfr = ecm_test(Y, X, COVID)
assert abs(F_ref - F_chk) < 1e-8 and abs(t_ref - t_chk) < 1e-8, (F_ref, F_chk, t_ref, t_chk)
NOBS = int(res.nobs)
cv = bounds_cv(NOBS, len(REGS), 3, nsim=40000, seed=42)
F_UP5, F_UP10, F_LO5 = cv["F"]["5"][1], cv["F"]["10"][1], cv["F"]["5"][0]
T_UP5 = cv["t"]["5"][1]
SIGMA_ECM = float(np.sqrt(res.scale))
e_ecm = res.resid
RHO_ECM = float(np.corrcoef(e_ecm[1:], e_ecm[:-1])[0, 1])

OUTD = {"design": dict(dep=DEP, regs=REGS, p=P, q=Q, nobs=NOBS, dfr=dfr, seed=SEED,
                       F_upper_5=F_UP5, F_upper_10=F_UP10, F_lower_5=F_LO5, t_upper_5=T_UP5,
                       check_F=F_chk, check_t=t_chk)}

# ------------------------------------------------------------ calibration
dy = np.diff(Y)
# null model: AR(1) in first differences with intercept and COVID dummy
Zn = np.column_stack([np.ones(T - 2), dy[:-1], COVID[2:]])
fit_n = sm.OLS(dy[1:], Zn).fit()
MU, GAM, DCOV = [float(v) for v in fit_n.params]
SIG0 = float(np.sqrt(fit_n.scale))
RES0 = fit_n.resid - fit_n.resid.mean()
# null model with regime-specific drift
Zr = np.column_stack([1 - REGIME[2:], REGIME[2:], dy[:-1], COVID[2:]])
fit_r = sm.OLS(dy[1:], Zr).fit()
MU1, MU2, GAMR, DCOVR = [float(v) for v in fit_r.params]
SIGR = float(np.sqrt(fit_r.scale))
# regressor processes
xpar = {}
innov = np.empty((T - 2, len(REGS)))
for j in range(len(REGS)):
    if j in I1:
        dx = np.diff(X[:, j])
        f = sm.OLS(dx[1:], sm.add_constant(dx[:-1])).fit()
    else:
        f = sm.OLS(X[2:, j], sm.add_constant(X[1:-1, j])).fit()
    xpar[j] = (float(f.params[0]), float(f.params[1]))
    innov[:, j] = f.resid
COVX = np.cov(innov, rowvar=False)
CHOLX = np.linalg.cholesky(COVX)

OUTD["calibration"] = dict(
    null_model=dict(mu=MU, gamma=GAM, covid=DCOV, sigma=SIG0),
    null_model_regime=dict(mu_writeoff=MU1, mu_expansion=MU2, gamma=GAMR, covid=DCOVR, sigma=SIGR),
    ecm_sigma=SIGMA_ECM, ecm_resid_ar1=RHO_ECM, ecm_bg_p=diagnostics(res)["lm"][1],
    regressors={REGS[j]: dict(kind="I(1)" if j in I1 else "I(0)", const=xpar[j][0], ar=xpar[j][1],
                              innov_sd=float(np.sqrt(COVX[j, j]))) for j in range(len(REGS))})
print("calibration:", json.dumps(OUTD["calibration"], indent=1)[:900])


def gen_y_null(errors, mu1=MU, mu2=MU, gam=GAM, dcov=DCOV):
    y = np.empty(T)
    y[0] = Y[0]
    y[1] = Y[1]
    for t in range(2, T):
        mu = mu2 if REGIME[t] == 1 else mu1
        y[t] = y[t - 1] + mu + gam * (y[t - 1] - y[t - 2]) + dcov * COVID[t] + errors[t]
    return y


def ar_errors(rho, sigma):
    """AR(1) errors with unconditional standard deviation sigma."""
    eta = rng.standard_normal(T) * sigma * np.sqrt(1 - rho ** 2)
    e = np.empty(T)
    e[0] = rng.standard_normal() * sigma
    for t in range(1, T):
        e[t] = rho * e[t - 1] + eta[t]
    return e


def gen_x_sim():
    v = rng.standard_normal((T, len(REGS))) @ CHOLX.T
    Xs = np.empty((T, len(REGS)))
    Xs[:2] = X[:2]
    for j in range(len(REGS)):
        a, b = xpar[j]
        for t in range(2, T):
            if j in I1:
                Xs[t, j] = Xs[t - 1, j] + a + b * (Xs[t - 1, j] - Xs[t - 2, j]) + v[t, j]
            else:
                Xs[t, j] = a + b * Xs[t - 1, j] + v[t, j]
    return Xs


def run_size(label, make, nrep=NREP_SIZE):
    F = np.empty(nrep); t = np.empty(nrep)
    for r in range(nrep):
        y, Xm = make()
        F[r], t[r], _ = ecm_test(y, Xm, COVID)
    out = dict(design=label, nrep=nrep, mean_F=float(F.mean()),
               rej_F_upper_5=float((F > F_UP5).mean()), rej_F_upper_10=float((F > F_UP10).mean()),
               rej_F_lower_5=float((F > F_LO5).mean()),
               rej_joint_5=float(((F > F_UP5) & (t < T_UP5)).mean()))
    print(f"{label:62} mean F={out['mean_F']:.2f}  F>U5: {out['rej_F_upper_5']:.3f}  "
          f"F>U10: {out['rej_F_upper_10']:.3f}  joint: {out['rej_joint_5']:.3f}")
    return out


print("\nSIZE (no level relationship), %d replications, nominal 5%%" % NREP_SIZE)
size = []

def ideal():
    e = rng.standard_normal((T, 5))
    Xi = np.empty((T, 4))
    Xi[:, :2] = np.cumsum(e[:, 1:3], axis=0)
    Xi[:, 2:] = e[:, 3:]
    return np.cumsum(e[:, 0]), Xi

size.append(run_size("A. Idealised: unit random walks, i.i.d. errors (revision 1 design)", ideal))
size.append(run_size("B. Calibrated null, observed regressors, Gaussian errors",
                     lambda: (gen_y_null(rng.standard_normal(T) * SIG0), X)))
size.append(run_size("C. Calibrated null, observed regressors, resampled residuals",
                     lambda: (gen_y_null(rng.choice(RES0, size=T, replace=True)), X)))
size.append(run_size("D. As B with AR(1) errors, rho = 0.3", lambda: (gen_y_null(ar_errors(0.3, SIG0)), X)))
size.append(run_size("E. As B with AR(1) errors, rho = 0.5", lambda: (gen_y_null(ar_errors(0.5, SIG0)), X)))
size.append(run_size("F. As B with a break in drift at 2022Q1 (estimated)",
                     lambda: (gen_y_null(rng.standard_normal(T) * SIGR, MU1, MU2, GAMR, DCOVR), X)))
size.append(run_size("G. As F with AR(1) errors, rho = 0.3",
                     lambda: (gen_y_null(ar_errors(0.3, SIGR), MU1, MU2, GAMR, DCOVR), X)))
size.append(run_size("H. Calibrated null, regressors simulated from estimated processes",
                     lambda: (gen_y_null(rng.standard_normal(T) * SIG0), gen_x_sim())))
OUTD["size"] = size


# ------------------------------------------------------------ power
def gen_y_alt(phi, beta, errors):
    """ECM towards y = c + beta * mortgage share, with the estimated COVID effect."""
    c = Y[0] - beta * X[0, 0]
    y = np.empty(T)
    y[0] = Y[0]
    for t in range(1, T):
        y[t] = y[t - 1] + phi * (y[t - 1] - c - beta * X[t - 1, 0]) + DCOV * COVID[t] + errors[t]
    return y


def run_power(phi, beta, rho=0.0, nrep=NREP_POWER):
    rF = 0; rJ = 0
    for r in range(nrep):
        e = ar_errors(rho, SIGMA_ECM) if rho else rng.standard_normal(T) * SIGMA_ECM
        F, t, _ = ecm_test(gen_y_alt(phi, beta, e), X, COVID)
        rF += F > F_UP5
        rJ += (F > F_UP5) and (t < T_UP5)
    return rF / nrep, rJ / nrep


print("\nPOWER (level relationship present), %d replications per cell, 5%% level" % NREP_POWER)
PHIS = (-0.10, -0.20, -0.30, -0.50, -0.75)
BETAS = (0.05, 0.10, 0.20, 0.40)
npl_end, npl_mean = float(d["npl"].iloc[-1]), float(d["npl"].mean())
power = {"phis": list(PHIS), "betas": list(BETAS), "sigma": SIGMA_ECM, "cells": [],
         "beta_in_npl_pp_end": {str(b): b * npl_end * (1 - npl_end / 100) for b in BETAS},
         "beta_in_npl_pp_mean": {str(b): b * npl_mean * (1 - npl_mean / 100) for b in BETAS},
         "signal_to_noise": {str(b): float(b * np.std(np.diff(X[:, 0]), ddof=1) / SIGMA_ECM) for b in BETAS}}
print(f"{'|phi|':>6}{'half-life':>10}" + "".join(f"{'beta='+str(b):>14}" for b in BETAS) + "   (F rule / joint F and t rule)")
for phi in PHIS:
    hl = float(np.log(0.5) / np.log(1 + phi))
    line = f"{abs(phi):6.2f}{hl:10.1f}"
    for b in BETAS:
        pF, pJ = run_power(phi, b)
        power["cells"].append(dict(phi=phi, half_life=hl, beta=b, power_F=float(pF), power_joint=float(pJ)))
        line += f"   {pF:.3f}/{pJ:.3f}"
    print(line)
# robustness: serially correlated equilibrium errors
power["ar_errors"] = []
for phi in PHIS:
    pF, pJ = run_power(phi, 0.10, rho=0.3)
    power["ar_errors"].append(dict(phi=phi, beta=0.10, rho=0.3, power_F=float(pF), power_joint=float(pJ)))
    print(f"  rho=0.3, beta=0.10, |phi|={abs(phi):.2f}: {pF:.3f}/{pJ:.3f}")
OUTD["power"] = power

jdump(OUTD, "size_power.json")
print("\nsaved -> output/size_power.json")
