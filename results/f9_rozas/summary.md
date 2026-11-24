# F9: escenario único de Homem-de-Mello et al. (2026), propuesto por H. Rozas

Fase A: 4050 consultas analizadas (1050 con P2 y P3) en 50 realizaciones (familia DET acelerada, 40 kg/s, N2 + pérdida esperada ρ = 1000, R = 0.1, fallas del DET).

## E0–E1. Reproducción de la consulta y supuestos

| Comprobación | Mediana | Peor caso |
|---|---|---|
| E0: |u_QP capturado − u_lazo|∞ [kg/s] | 4.6e-06 | 0.0113 |
| λ_min(GᵀG) (Asm. 3 sin R) | 1.6e-06 | 1.3e-06 |
| λ_min(GᵀG + R DᵀD) (Asm. 3 con R = 0.1) | 0.0022 | 0.0021 |
| máx. G_ki (k ≥ i) [K/(kg/s)] (Lema 1 iii: ≤ 0) | -0.0037 | -0.0034 |

Columnas de G con signo constante (orthant del Lema 1 iii): 100.0 % de las consultas.

## E2. Teorema 1, constructivo: el escenario ξ* reproduce el plan SAA

| Formulación | Caso del teorema | |u_ξ* − u_SAA|∞ mediana | p99 | máx. | ΔJ máx. |
|---|---|---|---|---|---|
| P1 pérdida esperada, ρ = 100 | (i) T fijo | 7.9e-06 | 0.00034 | 0.00088 | 3.6e-05 |
| P1 pérdida esperada, ρ = 1000 | (i) T fijo | 1.4e-05 | 9.3e-05 | 0.00036 | 2.4e-05 |
| P2 fallas parciales, ρ = 100, una fila (implementada) | fuera del Lema 1 (m=1<n) | 1.1e-05 | 0.0053 | 0.0076 | 3.2e-05 |
| P2 fallas parciales, ρ = 1000, una fila | fuera del Lema 1 (m=1<n) | 1.8e-05 | 0.0104 | 0.0168 | 0.00024 |
| P2 por paso (m = n), ρ = 100 | (ii) Lema 1 (iii) | 6.9e-05 | 0.0182 | 0.0255 | 9.5e-05 |
| P2 por paso (m = n), ρ = 1000 | (ii) Lema 1 (iii) | 1.7e-05 | 0.0163 | 0.0222 | 8.1e-05 |
| P3 pronóstico de fuga, L1 | (i) T fijo | 5.2e-05 | 0.00031 | 0.00051 | 3.2e-07 |

El ξ* de P1 es un margen distinto en cada paso (dispersión entre pasos: mediana 2.4 K): no es la curva de calentamiento de ninguna duración única.
El ξ* de P3, traducido a error de pronóstico de fuga, está a 1.8e+05 σ (mediana) del pronóstico: existe, pero fuera de toda realización plausible (el artículo lo anticipa: ξ* puede estar fuera del soporte).

## E3. Bilevel (38)–(39) sobre escenarios FÍSICOS (sin contexto, como §7.1 del artículo)

Pérdida = J_E(u_escenario) − J_E(u_SAA), evaluada con 1000 cuantiles de D; VSS = J_E(N2 sin riesgo) − J_E(u_SAA).

| Escenario único | ρ | |Δu|∞ mediana | |Δu₁| (acción aplicada) p99 | pérdida / VSS (agregada) | consultas con |Δu|∞ < 0.01 |
|---|---|---|---|---|---|
| duración óptima D̂* (bilevel) | 100 | 4.5e-06 | 2.38 | 0.0632 | 79 % |
| cuantil 0.98 de D | 100 | 2.4e-05 | — | 35.3 | 60 % |
| mediana de D | 100 | 1.8e-05 | — | 1 | 76 % |
| media de D (pronóstico puntual) | 100 | 1.7e-05 | — | 1 | 76 % |
| duración óptima D̂* (bilevel) | 1000 | 0.281 | 3.5 | 0.005 | 44 % |
| cuantil 0.98 de D | 1000 | 0.612 | — | 0.132 | 42 % |
| mediana de D | 1000 | 11.2 | — | 1 | 23 % |
| media de D (pronóstico puntual) | 1000 | 11.2 | — | 1 | 23 % |

ρ = 100: D̂* mediana 170 s (p10–p90 -1–196 s); margen δ̂* mediana 44.7 K.

ρ = 1000: D̂* mediana 187 s (p10–p90 -1–288 s); margen δ̂* mediana 48.3 K.

| P2 / P3: escenario físico | ρ | |Δu|∞ mediana | máx. | pérdida / VSS |
|---|---|---|---|---|
| P2: una falla (ĵ, η̂, D̂) | 100 | 1.2e-05 | 0.287 | 0.578 (VSS medio 0.00018) |
| P2: una falla (ĵ, η̂, D̂) | 1000 | 1.8e-05 | 2.87 | 0.637 (VSS medio 0.0182) |
| P3: equivalente cierto (δL = 0) | — | 1.6e-08 | 0.189 | ΔJ medio 7.2e-05 (J medio 205) |
| P3: δL = s·σ (bilevel en s) | — | 1.6e-08 | 0.192 | ΔJ medio 3.2e-05 |

Efecto del error de pronóstico sobre la PCT predicha: σ_y ≤ 0.0042 K.

## CPU frente a N (tabla 3 del artículo)

| N | SAA (un QP) [s] | búsqueda del escenario [s] | |Δu|∞ |
|---|---|---|---|
| 10 | 0.0237 | 2.67 | 6.6e-05 |
| 40 | 0.0578 | 2.64 | 2.4e-06 |
| 100 | 0.125 | 2.64 | 1.1e-06 |
| 400 | 0.514 | 2.65 | 2.5e-08 |
| 1000 | 2.11 | 2.66 | 1.9e-08 |
| 4000 | 8.55 | 2.73 | 1.8e-08 |

## E4. Caso contextual: pronóstico Φ(θ, x) del escenario (algoritmo 1 + Nelder–Mead)

Entrenamiento: semillas 201–235; prueba: 236–250. Costo esperado exacto (1000 cuantiles) por consulta de prueba.

| ρ | Política | arrepentimiento medio vs SAA | máx. | fracción del VSS capturada | |Δu₁| medio [kg/s] |
|---|---|---|---|---|---|
| 100 | nominal (risk ignored) | 7.64 | 222 | 0 | 0.21 |
| 100 | mean duration E[D] (point forecast) | 7.64 | 222 | 9.9e-07 | 0.21 |
| 100 | median duration | 7.64 | 222 | 1e-06 | 0.21 |
| 100 | constant margin (no features) | 1.46 | 57.5 | 0.809 | 0.0749 |
| 100 | constant duration D̂ (no features) | 5.29 | 146 | 0.308 | 0.17 |
| 100 | AD linear Φ(θ;x) → margin | 0.822 | 34.1 | 0.892 | 0.0465 |
| 100 | AD linear Φ(θ;x) → duration | 1.09 | 37.2 | 0.857 | 0.0528 |
| 100 | AD quadratic Φ(θ;x) → margin | 0.942 | 39.3 | 0.877 | 0.0515 |
| 100 | [E-loss] constant duration D̂ | 4.94 | 122 | 0.353 | 0.163 |
| 100 | [E-loss] linear Φ(θ;x) → margin | 0.634 | 9.82 | 0.917 | 0.0487 |
| 100 | [E-loss] linear Φ(θ;x) → duration | 0.856 | 28.1 | 0.888 | 0.0528 |
| 100 | oracle ξ*_x (best one scenario per query) | 0.552 | 9.5 | 0.928 | 0.0621 |
| 100 | SAA pi_S (200 quantiles; target) | 0 | 0 | 1 | 0 |
| 1000 | nominal (risk ignored) | 377 | 4.44e+03 | 0 | 2.42 |
| 1000 | mean duration E[D] (point forecast) | 377 | 4.44e+03 | -4.6e-07 | 2.42 |
| 1000 | median duration | 377 | 4.44e+03 | -1e-08 | 2.42 |
| 1000 | constant margin (no features) | 51.1 | 1.01e+03 | 0.864 | 1.57 |
| 1000 | constant duration D̂ (no features) | 40.7 | 1.14e+03 | 0.892 | 1.3 |
| 1000 | AD linear Φ(θ;x) → margin | 11.6 | 483 | 0.969 | 1.09 |
| 1000 | AD linear Φ(θ;x) → duration | 32.2 | 730 | 0.915 | 1.45 |
| 1000 | AD quadratic Φ(θ;x) → margin | 16.7 | 495 | 0.956 | 1.07 |
| 1000 | [E-loss] constant duration D̂ | 40 | 824 | 0.894 | 1.23 |
| 1000 | [E-loss] linear Φ(θ;x) → margin | 1.27 | 20.7 | 0.997 | 0.989 |
| 1000 | [E-loss] linear Φ(θ;x) → duration | 1.55 | 36.4 | 0.996 | 0.982 |
| 1000 | oracle ξ*_x (best one scenario per query) | 0.207 | 18.5 | 0.999 | 0.917 |
| 1000 | SAA pi_S (200 quantiles; target) | 0 | 0 | 1 | 0 |

Curva de aprendizaje (teorema 2: θ*_N converge al crecer N):

| ρ | N entrenamiento | fracción del VSS capturada (media de réplicas) |
|---|---|---|
| 100 | 25 | -0.8 |
| 100 | 100 | 0.435 |
| 100 | 300 | 0.469 |
| 100 | 1000 | 0.773 |
| 100 | 2835 | 0.892 |
| 1000 | 25 | 0.553 |
| 1000 | 100 | 0.512 |
| 1000 | 300 | 0.761 |
| 1000 | 1000 | 0.929 |
| 1000 | 2835 | 0.969 |

## E5. Lazo cerrado (sin fallas, Rao–Blackwell μ×10⁴, 40 kg/s, R = 0.1)

| Caso | n | P(1220 K | falla) | pre-enfriamiento medio [K] | masa [kg] | variación total [kg/s] | t p99 [s] |
|---|---|---|---|---|---|---|---|
| N2 | 50 | 0.0311 | -0.0273 | 3.72e+03 | 60.7 | 0.0485 |
| SAA_r100 | 50 | 0.0258 | 1.8 | 3.86e+03 | 66.9 | 0.293 |
| SAA_r1000 | 50 | 0.0137 | 9.41 | 4.51e+03 | 69.7 | 0.299 |
| scen_ADlinDE_r100 | 50 | 0.0198 | 4.7 | 4.1e+03 | 2.24e+03 | 0.0324 |
| scen_ADlinDE_r1000 | 50 | 0.0145 | 8.43 | 4.42e+03 | 1.54e+03 | 0.0363 |
| scen_ADlinE_r100 | 50 | 0.0259 | 1.73 | 3.85e+03 | 74.1 | 0.0306 |
| scen_ADlinE_r1000 | 50 | 0.0147 | 8.22 | 4.4e+03 | 1.69e+03 | 0.0321 |
| scen_ADlin_r100 | 50 | 0.024 | 2.65 | 3.92e+03 | 150 | 0.0501 |
| scen_ADlin_r1000 | 50 | 0.0148 | 8.19 | 4.4e+03 | 1.69e+03 | 0.0475 |
| scen_ADslowE_r100 | 50 | 0.0257 | 1.83 | 3.86e+03 | 71.9 | 0.0308 |
| scen_ADslowE_r1000 | 50 | 0.0148 | 8.19 | 4.4e+03 | 84.3 | 0.0304 |
| scen_Dconst_r100 | 50 | 0.0294 | 0.364 | 3.75e+03 | 63.6 | 0.0322 |
| scen_Dconst_r1000 | 50 | 0.0178 | 5.6 | 4.17e+03 | 84 | 0.0537 |

