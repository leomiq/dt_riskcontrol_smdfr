"""Summary tables of the F6 paper campaign (plan §10.1; julia/f6_campaign.jl).

Reads results/f6_paper/<campaign>.csv (and *_audit.csv) and writes results/f6_paper/summary.md:
  sens  P(cross | failure) [Rao–Blackwell] by ACS capacity and failure acceleration μ×, for
        N2, PI, PID and the risk bounds designed with the matched μ; robustness of the bounds
        designed with μ×1e4 when the true μ differs; pre-cooling (K below y_ref) and mass
  form  every risk formulation on the risk–pre-cooling frontier (40 kg/s, μ×1e4)
  eta   P(cross | failure) for the four severity-probability vectors
  det, nom, ext  thesis control metrics per controller, crossings of 1220 K, solve times,
        and the per-query optimality audit (gap of the loop decision vs the best N3)
Paired comparisons use common random numbers (same seed); intervals are mean ± 1.96 s.e.
Only csv + numpy (pandas is blocked on leo_acer).

Usage: python python/f6_summary.py [results_dir=results/f6_paper]
"""
import csv, sys, os, math
from collections import defaultdict
import numpy as np

D = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), "..", "results", "f6_paper")
out = []
say = out.append


def load(name):
    p = os.path.join(D, name + ".csv")
    if not os.path.isfile(p):
        return []
    with open(p, newline="") as f:
        return list(csv.DictReader(f))


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return math.nan


def ci(v):
    v = np.asarray([x for x in v if np.isfinite(x)], float)
    if len(v) == 0:
        return "—"
    if len(v) == 1:
        return f"{v[0]:.4g}"
    return f"{v.mean():.4g} ± {1.96 * v.std(ddof=1) / math.sqrt(len(v)):.2g}"


def med(v):
    v = np.asarray([x for x in v if np.isfinite(x)], float)
    return f"{np.median(v):.3g}" if len(v) else "—"


def table(head, rows):
    say("| " + " | ".join(head) + " |")
    say("|" + "---|" * len(head))
    for r in rows:
        say("| " + " | ".join(str(c) for c in r) + " |")
    say("")


def by(rows, *keys):
    g = defaultdict(list)
    for r in rows:
        g[tuple(r[k] for k in keys)].append(r)
    return g


def col(rows, k):
    return [fnum(r[k]) for r in rows]


def paired(rows, ref_rows, k):
    ref = {r["seed"]: fnum(r[k]) for r in ref_rows}
    return [fnum(r[k]) - ref[r["seed"]] for r in rows if r["seed"] in ref]


# ------------------------------------------------------------------ sens
S = load("sens")
if S:
    say("## sens: failure acceleration μ× and ACS capacity (Rao–Blackwell, P(cross | failure))\n")
    say(f"{len({r['seed'] for r in S})} realisations. Rows: controller (and the μ× it was designed with); "
        "columns: true μ×. Matched design in bold. Pre-cooling = K below y_ref without failures.\n")
    G = by(S, "act", "case")
    for act in sorted({r["act"] for r in S}, key=lambda a: float(a[3:])):
        say(f"### {act}\n")
        rows = []
        base = G.get((act, "N2"), [])
        for case in sorted({r["case"] for r in S if r["act"] == act}):
            R = G[(act, case)]
            mu_d = case.split("_mu")[-1] if "_mu" in case else ""
            cells = []
            for m in ("3e3", "1e4", "3e4"):
                c = ci(col(R, f"p_cond_mu{m}"))
                cells.append(f"**{c}**" if m == mu_d else c)
            pre = ci([-x for x in col(R, "mean_dev")])
            mass = np.nanmean(col(R, "mass_plant")) / np.nanmean(col(base, "mass_plant")) - 1 if base else math.nan
            rows.append([case] + cells + [pre, f"{100 * mass:+.0f} %", ci(col(R, "IAE"))])
        table(["case", "μ×3e3", "μ×1e4", "μ×3e4", "pre-cooling [K]", "mass vs N2", "IAE [K·s]"], rows)

# ------------------------------------------------------------------ form
F = load("form") + [r for r in S if r["act"] == "cap40" and (r["case"] in ("N2", "PI", "PID") or r["case"].endswith("_mu1e4"))]
if F:
    say("## form: risk–pre-cooling frontier (40 kg/s, μ×1e4)\n")
    G = by(F, "case"); base = G.get(("N2",), [])
    rows = []
    for (case,), R in sorted(G.items(), key=lambda kv: -np.nanmean(col(kv[1], "p_cond_mu1e4"))):
        rows.append([case, len(R), ci(col(R, "p_cond_mu1e4")), ci([-x for x in col(R, "mean_dev")]),
                     f"{100 * (np.nanmean(col(R, 'mass_plant')) / np.nanmean(col(base, 'mass_plant')) - 1):+.0f} %" if base else "—",
                     med(col(R, "t_p99"))])
    table(["case", "n", "P(cross | failure)", "pre-cooling [K]", "mass vs N2", "p99 solve [s]"], rows)

# ------------------------------------------------------------------ eta
E = load("eta")
if E:
    say("## eta: partial failures, severity probabilities p_η over η = (0.75, 0.5, 0.25, 0)\n")
    say("base (.4,.3,.2,.1), uniform (.25 each), pess (.1,.2,.3,.4), twotrain (0,.7,0,.3). "
        "`total only` = total failures only (p_0 = 1).\n")
    rows = []
    for (case,), R in sorted(by(E, "case").items()):
        rows.append([case, len(R)] + [ci(col(R, f"p_cond_{p}")) for p in ("base", "uniform", "pess", "twotrain")]
                    + [ci(col(R, "p_total_only")), ci([-x for x in col(R, "mean_dev")])])
    table(["case", "n", "base", "uniform", "pess", "twotrain", "total only", "pre-cooling [K]"], rows)

# ------------------------------------------------------------------ closed-loop families
TH = [("IAE", "IAE [K·s]"), ("settling_time", "settling [s]"), ("min_overshoot", "min overshoot [K]"),
      ("steady_state_error", "ss error [K]"), ("total_variation", "TV [kg/s]"), ("min_safety_margin", "min DSM (1244) [K]"),
      ("peak_PCT", "peak PCT [K]")]
for camp, title in (("det", "closed loop with DET-accelerated failures (40 kg/s)"),
                    ("nom", "nominal family (thesis setting)"), ("ext", "extreme family (F5a.3)")):
    C = load(camp)
    if not C:
        continue
    say(f"## {camp}: {title}\n")
    say("Thesis metrics (Table 10 of the control paper; `thesis_metrics`): mean ± 95 % CI; settling = "
        "median (Inf = never settles in the 0.30 K band), with the fraction that settles. "
        "Crossings = runs that reach 1220 K.\n")
    for act in sorted({r["act"] for r in C}):
        Ca = [r for r in C if r["act"] == act]
        say(f"### actuator {act} ({len({r['seed'] for r in Ca})} realisations)\n")
        G = by(Ca, "case"); ref = G.get(("N2b",), [])
        rows = []
        for (case,), R in sorted(G.items()):
            st = col(R, "settling_time")
            row = [case, len(R), f"{sum(r['limit'] == 'true' for r in R)}/{len(R)}"]
            for k, _ in TH:
                row.append(f"{med(st)} ({np.mean(np.isfinite(st)):.0%})" if k == "settling_time" else ci(col(R, k)))
            row += [ci(paired(R, ref, "IAE")) if ref and case != "N2b" else "—", med(col(R, "t_p99"))]
            rows.append(row)
        table(["case", "n", "crossings"] + [h for _, h in TH] + ["ΔIAE vs N2b (paired)", "p99 solve [s]"], rows)
    A = load(camp + "_audit")
    if A:
        say("**Per-query optimality audit** (P0 cost J on the ode1be ROM, horizon K = 26; J_best = best of the "
            "loop decision, N1, N2, N3 warm/polished/cold). In most queries the set-point is tracked exactly and "
            "J_best ≈ 1e-12, so the relative gap J/J_best − 1 is meaningless there: the excess is reported as an "
            "RMS tracking error over the horizon, ΔRMS = √(J/K) − √(J_best/K) [K], and the relative gap only "
            "where J_best ≥ 1 K².\n")
        rows = []
        for (key_case,), R in sorted(by([dict(r, kc=r["key"].split("|")[3] + "|" + r["key"].split("|")[4]) for r in A], "kc").items()):
            Jb = np.array(col(R, "Jbest")); g = np.array(col(R, "gap_loop")); gc = np.array(col(R, "gap_N3cold"))
            drms = np.sqrt(Jb * (1 + g) / 26) - np.sqrt(Jb / 26)
            dcold = np.sqrt(Jb * (1 + gc) / 26) - np.sqrt(Jb / 26)
            big = Jb >= 1
            rows.append([key_case, len(R), f"{np.mean(Jb < 1e-2):.0%}", f"{np.quantile(drms, 0.99):.1e}", f"{drms.max():.1e}",
                         f"{np.median(g[big]):.1e} / {g[big].max():.1e} (n={big.sum()})" if big.any() else "—",
                         f"{np.mean(dcold > 1e-3):.1%}"])
        table(["actuator|loop", "audited queries", "J_best < 0.01 K²", "p99 ΔRMS [K]", "max ΔRMS [K]",
               "rel. gap where J_best ≥ 1: median / max", "N3 cold start worse by > 1 mK RMS"], rows)

for camp in ("sens", "form", "eta", "det", "nom", "ext"):
    p = os.path.join(D, camp + "_errors.log")
    if os.path.isfile(p) and os.path.getsize(p) > 0:
        say(f"**{camp}: errors logged in {camp}_errors.log** ({open(p).read().count('[20')} tasks)\n")

dst = os.path.join(D, "summary.md")
with open(dst, "w", encoding="utf-8") as f:
    f.write("# F6 paper campaign — summary\n\n" + "\n".join(out))
print("written", dst)
