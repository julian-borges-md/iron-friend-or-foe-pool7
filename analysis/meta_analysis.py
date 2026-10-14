"""
Iron Overload and Male Hypogonadism: Random-Effects Meta-Analysis (Python cross-check)
TDT male hypogonadism prevalence. PROSPERO CRD420261344877. Borges JYV.

Model: DerSimonian-Laird. Transformation: Freeman-Tukey double arcsine.
Back-transform: Miller (1978), harmonic mean n.
Dependencies: numpy, scipy, pandas.
"""
import numpy as np
import pandas as pd
from scipy import stats

dat = pd.read_csv("data/tdt_male_hypogonadism_prevalence.csv")
n = dat["n_males"].to_numpy(float)
x = dat["events"].to_numpy(float)

y = 0.5 * (np.arcsin(np.sqrt(x / (n + 1))) + np.arcsin(np.sqrt((x + 1) / (n + 1))))
v = 1.0 / (4 * n + 2)
w = 1.0 / v

y_fe = np.sum(w * y) / np.sum(w)
Q = np.sum(w * (y - y_fe) ** 2)
df = len(n) - 1
C = np.sum(w) - np.sum(w ** 2) / np.sum(w)
tau2 = max(0.0, (Q - df) / C)
I2 = max(0.0, (Q - df) / Q) * 100

w_re = 1.0 / (v + tau2)
y_re = np.sum(w_re * y) / np.sum(w_re)
se_re = np.sqrt(1.0 / np.sum(w_re))
se_fe = np.sqrt(1.0 / np.sum(w))

n_h = len(n) / np.sum(1.0 / n)


def back_transform(t: float) -> float:
    s = np.sin(2 * t)
    return 0.5 * (1 - np.sign(np.cos(2 * t)) * np.sqrt(1 - (s + (s - 1 / s) / n_h) ** 2))


print(f"Random effects (DL): {back_transform(y_re)*100:.1f}% "
      f"[{back_transform(y_re-1.96*se_re)*100:.1f} to {back_transform(y_re+1.96*se_re)*100:.1f}]")
print(f"Fixed effects (IV):  {back_transform(y_fe)*100:.1f}% "
      f"[{back_transform(y_fe-1.96*se_fe)*100:.1f} to {back_transform(y_fe+1.96*se_fe)*100:.1f}]")
print(f"Q = {Q:.2f}, df = {df}, p = {1-stats.chi2.cdf(Q, df):.2e}, I2 = {I2:.1f}%, tau2 = {tau2:.4f}")
print("RE weights (%):", dict(zip(dat["study"], np.round(w_re / np.sum(w_re) * 100, 1))))
