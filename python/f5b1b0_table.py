"""LaTeX table of the thesis indices for step 3.i.0 (docs/estado/11g_dsm_lambda_tabla.tex) from
results/f5b1b0_dsm_lambda/runs.csv. Means over realisations; settling time over the runs that
settle, with the count of runs that never settle in the 0.30 K band.
Usage: .venv/Scripts/python.exe python/f5b1b0_table.py
"""
import csv
import os
from collections import defaultdict
import numpy as np

ROOT = os.path.join(os.path.dirname(__file__), "..")
rows = list(csv.DictReader(open(os.path.join(ROOT, "results", "f5b1b0_dsm_lambda", "runs.csv"))))
BLOCKS = [("archived", "arch", "61", "Realización archivada, grilla de la tesis (61 consultas), bomba archivada 12.9~kg/s"),
          ("archived", "wide", "61", "Realización archivada, grilla de la tesis (61 consultas), bomba ampliada 43.6~kg/s"),
          ("nominal", "arch", "121", "Familia nominal, 20 realizaciones, bomba archivada"),
          ("nominal", "wide", "121", "Familia nominal, 20 realizaciones, bomba ampliada"),
          ("detacc", "cap40", "481", "Familia DET acelerada (con fallas), 40~kg/s")]
CASES = [("N2", "N2+sesgo"), ("N2_DSM_U1244_l10000", r"N2+DSM, $\lambda$=$10^4$, $U$=1244"),
         ("N2_DSM_U1220_l1000000", r"N2+DSM, $\lambda$=$10^6$, $U$=1220"), ("PID", "PID de la tesis"),
         ("PID_free", "PID de la tesis, sin límite"), ("PI", "PI del DET")]
out = [r"\begin{table}[H]\centering\scriptsize\setlength{\tabcolsep}{3pt}",
       r"\begin{tabularx}{\linewidth}{@{}Y r r r r r r r r@{}}\toprule",
       r"Controlador & IAE [K$\cdot$s] & Alcanza [s] & Asienta [s] & Err.\ estac.\ [K] & Sobrepaso mín.\ [K] & PCT pico [K] & VT [kg/s] & Llegan a 1220\\ \midrule"]
for st, act, n, title in BLOCKS:
    by = defaultdict(dict)
    for r in rows:
        if r["set"] == st and r["act"] == act and r["n_req"] == n:
            by[r["case"]][r["seed"]] = r
    if not by:
        continue
    nreal = len(by.get("N2", {}))
    out.append(r"\multicolumn{9}{@{}l}{\textit{" + title + (f" ({nreal} realizaciones)" if st == "detacc" else "") + r"}}\\")
    for c, lab in CASES:
        if c not in by:
            continue
        S = by[c]; v = lambda k: np.array([float(S[s][k]) for s in S])
        reach, sett = v("reach_time"), v("settling_time")
        fr, fs = reach[np.isfinite(reach)], sett[np.isfinite(sett)]
        nos = int(np.sum(~np.isfinite(sett)))
        rs = f"{fr.mean():.1f}" if fr.size else r"$\infty$"
        ss = (f"{fs.mean():.1f}" if fs.size else r"$\infty$") + (f" ({nos} no)" if nos and fs.size else "")
        lim = sum(S[s]["limit"] == "true" for s in S)
        out.append(f"{lab} & {v('IAE').mean():.1f} & {rs} & {ss} & {v('steady_state_error').mean():+.2f} & "
                   f"{v('min_overshoot').mean():+.2f} & {v('peak_PCT').mean():.2f} & {v('tv_plant').mean():.0f} & {lim}/{len(S)}\\\\")
    out.append(r"\addlinespace")
out[-1] = r"\bottomrule"
out += [r"\end{tabularx}",
        r"\caption{Métricas de la tesis (medias sobre las realizaciones). «Alcanza»: primer instante con "
        r"$|e|\le0.30$~K; «Asienta»: instante desde el cual no sale de la banda ($\infty$: nunca; entre "
        r"paréntesis, las realizaciones que no se asientan). Error estacionario: media de $y_{\rm ref}-\mathrm{PCT}$ "
        r"en el último 20\,\%. VT: variación total del caudal entregado. Tiempos desde la primera consulta. "
        r"Las demás combinaciones de $\lambda_{\rm dsm}$ y $U$ dan lo mismo que N2 (\code{results/f5b1b0\_dsm\_lambda/runs.csv}).}",
        r"\label{tab:dsm-lambda}", r"\end{table}"]
open(os.path.join(ROOT, "docs", "estado", "11g_dsm_lambda_tabla.tex"), "w", encoding="utf-8").write("\n".join(out) + "\n")
print("\n".join(out))
