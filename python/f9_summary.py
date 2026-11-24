"""Summary of F9 (results/f9_rozas/): the one-scenario result of Homem-de-Mello et al. (2026),
proposed by H. Rozas, checked on the SMDFR MPC queries with the paper's methodology.

Reads A_queries.csv (phase A, one row per analysed query), B_results.csv / B_learning.csv
(phase B, contextual forecasts) and C_rb.csv (phase C, closed loop) if present, and writes
results/f9_rozas/summary.md. Usage: python python/f9_summary.py
"""
import os, csv
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results", "f9_rozas")


def rows(name):
    p = os.path.join(RES, name + ".csv")
    return list(csv.DictReader(open(p, encoding="utf-8"))) if os.path.isfile(p) else []


def col(R, k):
    v = []
    for r in R:
        try:
            v.append(float(r.get(k, "nan")))
        except ValueError:
            v.append(np.nan)
    v = np.array(v)
    return v[np.isfinite(v)]


def q(v, p):
    return np.percentile(v, p) if len(v) else np.nan


def g(x):
    return "—" if not np.isfinite(x) else f"{x:.2g}" if abs(x) < 1e-2 or abs(x) >= 1e4 else f"{x:.3g}"


L = ["# F9: escenario único de Homem-de-Mello et al. (2026), propuesto por H. Rozas", ""]
A = rows("A_queries")
if A:
    F = [r for r in A if r["full"] == "true"]
    seeds = sorted({r["seed"] for r in A})
    L += [f"Fase A: {len(A)} consultas analizadas ({len(F)} con P2 y P3) en {len(seeds)} realizaciones "
          "(familia DET acelerada, 40 kg/s, N2 + pérdida esperada ρ = 1000, R = 0.1, fallas del DET).", ""]
    # ---------------- E0 / E1
    L += ["## E0–E1. Reproducción de la consulta y supuestos", "",
          "| Comprobación | Mediana | Peor caso |", "|---|---|---|"]
    for k, name in (("E0_gap", "E0: |u_QP capturado − u_lazo|∞ [kg/s]"),
                    ("lmin_GtG", "λ_min(GᵀG) (Asm. 3 sin R)"), ("lmin_H", "λ_min(GᵀG + R DᵀD) (Asm. 3 con R = 0.1)"),
                    ("G_maxpos", "máx. G_ki (k ≥ i) [K/(kg/s)] (Lema 1 iii: ≤ 0)")):
        v = col(A, k)
        worst = v.max() if k in ("E0_gap", "G_maxpos") else v.min()
        L.append(f"| {name} | {g(np.median(v))} | {g(worst)} |")
    ok = np.mean([r["G_colsign_ok"] == "true" for r in A])
    L += ["", f"Columnas de G con signo constante (orthant del Lema 1 iii): {100*ok:.1f} % de las consultas.", ""]
    # ---------------- E2
    L += ["## E2. Teorema 1, constructivo: el escenario ξ* reproduce el plan SAA", "",
          "| Formulación | Caso del teorema | |u_ξ* − u_SAA|∞ mediana | p99 | máx. | ΔJ máx. |", "|---|---|---|---|---|---|"]
    for k, name, case, R in (("r100_T1", "P1 pérdida esperada, ρ = 100", "(i) T fijo", A),
                             ("r1000_T1", "P1 pérdida esperada, ρ = 1000", "(i) T fijo", A),
                             ("p2r100_row1", "P2 fallas parciales, ρ = 100, una fila (implementada)", "fuera del Lema 1 (m=1<n)", F),
                             ("p2r1000_row1", "P2 fallas parciales, ρ = 1000, una fila", "fuera del Lema 1 (m=1<n)", F),
                             ("p2r100_lemma", "P2 por paso (m = n), ρ = 100", "(ii) Lema 1 (iii)", F),
                             ("p2r1000_lemma", "P2 por paso (m = n), ρ = 1000", "(ii) Lema 1 (iii)", F),
                             ("p3_T1", "P3 pronóstico de fuga, L1", "(i) T fijo", F)):
        v = col(R, k + "_gap"); d = col(R, k + "_dJ")
        L.append(f"| {name} | {case} | {g(np.median(v))} | {g(q(v, 99))} | {g(v.max())} | {g(np.abs(d).max())} |")
    sp = col(A, "r1000_T1_spread"); z = col(F, "p3_T1_zmax")
    L += ["", f"El ξ* de P1 es un margen distinto en cada paso (dispersión entre pasos: mediana {g(np.median(sp))} K): "
          "no es la curva de calentamiento de ninguna duración única.",
          f"El ξ* de P3, traducido a error de pronóstico de fuga, está a {g(np.median(z))} σ (mediana) del pronóstico: "
          "existe, pero fuera de toda realización plausible (el artículo lo anticipa: ξ* puede estar fuera del soporte).", ""]
    # ---------------- E3
    L += ["## E3. Bilevel (38)–(39) sobre escenarios FÍSICOS (sin contexto, como §7.1 del artículo)", "",
          "Pérdida = J_E(u_escenario) − J_E(u_SAA), evaluada con 1000 cuantiles de D; VSS = J_E(N2 sin riesgo) − J_E(u_SAA).", "",
          "| Escenario único | ρ | |Δu|∞ mediana | |Δu₁| (acción aplicada) p99 | pérdida / VSS (agregada) | consultas con |Δu|∞ < 0.01 |",
          "|---|---|---|---|---|---|"]
    for ρ in ("100", "1000"):
        vss = col(A, f"r{ρ}_Jnom").sum() - col(A, f"r{ρ}_JS").sum()
        for k, name in (("BL", "duración óptima D̂* (bilevel)"), ("q98", "cuantil 0.98 de D"),
                        ("median", "mediana de D"), ("mean", "media de D (pronóstico puntual)")):
            gap = col(A, f"r{ρ}_{k}_gap"); dJ = col(A, f"r{ρ}_{k}_dJ")
            g1 = col(A, f"r{ρ}_{k}_gap1") if k == "BL" else np.array([np.nan])
            L.append(f"| {name} | {ρ} | {g(np.median(gap))} | {g(q(g1, 99))} | {g(dJ.sum() / vss)} | "
                     f"{100*np.mean(gap < 0.01):.0f} % |")
    for ρ in ("100", "1000"):
        D = col(A, f"r{ρ}_BL_D"); dl = col(A, f"r{ρ}_BL_delta")
        L.append(f"\nρ = {ρ}: D̂* mediana {g(np.median(D[np.isfinite(D)]))} s (p10–p90 {g(q(D,10))}–{g(q(D,90))} s); "
                 f"margen δ̂* mediana {g(np.median(dl))} K.")
    L += ["", "| P2 / P3: escenario físico | ρ | |Δu|∞ mediana | máx. | pérdida / VSS |", "|---|---|---|---|---|"]
    for ρ in ("100", "1000"):
        vss = col(F, f"p2r{ρ}_VSS").sum()
        L.append(f"| P2: una falla (ĵ, η̂, D̂) | {ρ} | {g(np.median(col(F, f'p2r{ρ}_phys_gap')))} | "
                 f"{g(col(F, f'p2r{ρ}_phys_gap').max())} | {g(col(F, f'p2r{ρ}_phys_dJ').sum() / vss)} (VSS medio {g(vss/len(F))}) |")
    L.append(f"| P3: equivalente cierto (δL = 0) | — | {g(np.median(col(F, 'p3_CE_gap')))} | {g(col(F, 'p3_CE_gap').max())} | "
             f"ΔJ medio {g(np.mean(col(F, 'p3_CE_dJ')))} (J medio {g(np.mean(col(F, 'p3_J')))}) |")
    L.append(f"| P3: δL = s·σ (bilevel en s) | — | {g(np.median(col(F, 'p3_BL_gap')))} | {g(col(F, 'p3_BL_gap').max())} | "
             f"ΔJ medio {g(np.mean(col(F, 'p3_BL_dJ')))} |")
    L += ["", f"Efecto del error de pronóstico sobre la PCT predicha: σ_y ≤ {g(col(F, 'p3_sigma_y').max())} K.", ""]
    # ---------------- CPU
    S = [r for r in A if r.get("cpu_saa_N10", "nan") not in ("", "nan", "NaN")]
    if S:
        L += ["## CPU frente a N (tabla 3 del artículo)", "", "| N | SAA (un QP) [s] | búsqueda del escenario [s] | |Δu|∞ |",
              "|---|---|---|---|"]
        for N in (10, 40, 100, 400, 1000, 4000):
            L.append(f"| {N} | {g(np.median(col(S, f'cpu_saa_N{N}')))} | {g(np.median(col(S, f'cpu_bl_N{N}')))} | "
                     f"{g(np.median(col(S, f'cpu_gap_N{N}')))} |")
        L.append("")
B = rows("B_results")
if B:
    L += ["## E4. Caso contextual: pronóstico Φ(θ, x) del escenario (algoritmo 1 + Nelder–Mead)", "",
          "Entrenamiento: semillas 201–235; prueba: 236–250. Costo esperado exacto (1000 cuantiles) por consulta de prueba.", "",
          "| ρ | Política | arrepentimiento medio vs SAA | máx. | fracción del VSS capturada | |Δu₁| medio [kg/s] |",
          "|---|---|---|---|---|---|"]
    for r in B:
        L.append(f"| {r['rho']} | {r['policy']} | {g(float(r['regret_vs_SAA']))} | {g(float(r['regret_max']))} | "
                 f"{g(float(r['VSS_captured']))} | {g(float(r['u1_mae_vs_SAA']))} |")
    Lr = rows("B_learning")
    if Lr:
        L += ["", "Curva de aprendizaje (teorema 2: θ*_N converge al crecer N):", "",
              "| ρ | N entrenamiento | fracción del VSS capturada (media de réplicas) |", "|---|---|---|"]
        for ρ in sorted({r["rho"] for r in Lr}, key=float):
            for N in sorted({r["N_train"] for r in Lr if r["rho"] == ρ}, key=int):
                v = col([r for r in Lr if r["rho"] == ρ and r["N_train"] == N], "VSS_captured")
                L.append(f"| {ρ} | {N} | {g(v.mean())} |")
    L.append("")
C = rows("C_rb")
if C:
    L += ["## E5. Lazo cerrado (sin fallas, Rao–Blackwell μ×10⁴, 40 kg/s, R = 0.1)", "",
          "| Caso | n | P(1220 K | falla) | pre-enfriamiento medio [K] | masa [kg] | variación total [kg/s] | t p99 [s] |",
          "|---|---|---|---|---|---|---|---|"]
    cases = sorted({r["case"] for r in C}, key=lambda c: (c != "N2", c))
    for c in cases:
        R = [r for r in C if r["case"] == c]
        L.append(f"| {c} | {len(R)} | {g(col(R,'p_cond').mean())} | {g(-col(R,'mean_dev').mean())} | "
                 f"{g(col(R,'mass').mean())} | {g(col(R,'tv').mean())} | {g(q(col(R,'t_p99'),50))} |")
    L.append("")
open(os.path.join(RES, "summary.md"), "w", encoding="utf-8").write("\n".join(L) + "\n")
print("\n".join(L))
