"""F5a.2: paired summary of the closed-loop Monte Carlo (results/f5a2_monte_carlo/runs.csv).

Every controller faces the same realisations (common random numbers), so differences are
computed per realisation against a reference controller (N3 with the bias correction and
the same actuator configuration) and summarised with a t-based 95 % interval.

Usage: python python/f5a2_summary.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
df = pd.read_csv(ROOT / "results" / "f5a2_monte_carlo" / "runs.csv")
df["bias"] = df["bias"].astype(str).str.lower() == "true"
df["label"] = df["act"] + " " + df["ctrl"] + np.where(df["bias"], "+bias", "")

print("Realisations:")
real = df.drop_duplicates("seed")[["seed", "family", "A", "B", "drift", "t_act", "t_start", "eta",
                                    "fail_start", "fail_len", "fbias"]]
print(real.to_string(index=False, float_format=lambda v: f"{v:.3f}"))


def ci(x):
    x = np.asarray(x, float)
    if len(x) < 2:
        return np.nan
    return stats.t.ppf(0.975, len(x) - 1) * x.std(ddof=1) / np.sqrt(len(x))


rows = []
for fam, F in df.groupby("family"):
    for lab, S in F.groupby("label"):
        act = S["act"].iloc[0]
        ref_act = "phys" if act == "none" else act
        ref = F[(F["act"] == ref_act) & (F["ctrl"] == "N3") & F["bias"]].set_index("seed")
        s = S.set_index("seed")
        common = s.index.intersection(ref.index)
        d_ise = s.loc[common, "ise"] - ref.loc[common, "ise"]
        d_rel = d_ise / ref.loc[common, "ise"]
        rows.append(dict(family=fam, controller=lab, n=len(S),
                         ISE=S["ise"].mean(),
                         dISE=d_ise.mean(), dISE_ci=ci(d_ise), dISE_rel_pct=100 * d_rel.mean(),
                         PCT_max=S["pct_max"].mean(), t_above_s=S["t_above"].mean(),
                         err_last10=S["err_last10"].mean(), mass_plant=S["mass_plant"].mean(),
                         solve_p99_max=S["t_p99"].max(), fails=int(S["fails"].sum())))
out = pd.DataFrame(rows)
order = {"none zero": 0, "arch max": 1}
out["o"] = out["controller"].map(lambda c: order.get(c, 2 if c.startswith("arch") else 3))
out = out.sort_values(["family", "o", "controller"]).drop(columns="o")
pd.set_option("display.width", 200)
print("\nPaired against N3+bias with the same actuator (Δ = controller − reference, mean ± 95% CI):")
print(out.to_string(index=False, float_format=lambda v: f"{v:.4g}"))
out.to_csv(ROOT / "results" / "f5a2_monte_carlo" / "summary.csv", index=False)
print("\nwritten results/f5a2_monte_carlo/summary.csv")
