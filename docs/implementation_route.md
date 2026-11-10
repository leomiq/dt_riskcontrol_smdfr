# Plan de mejoras y vías de exploración: control predictivo informado por riesgo

Versión: 22 de octubre de 2026. Punto de partida: v1.0.0 (commit `3266149`).
Documento compañero de la presentación [`presentations/usm_risk_control.pdf`](../presentations/usm_risk_control.pdf);
usa la misma notación y la misma numeración de problemas (P0, P1a, P1b), vías (V1–V10) y etapas (E0–E5).
Actualizado tras la reunión técnica del 22-10-2026 con H. Rozas (U. de Chile) y R. Mena (USM).

---

## 1. En una página

**Qué existe.** Un controlador predictivo (MPC) que, cada 0.5 s durante una pérdida de
refrigerante (LOCA) en el SMDFR, elige la corriente de la bomba MHD del sistema de
enfriamiento auxiliar (ACS). Predice la temperatura pico de vaina (PCT) con un gemelo
digital *grey-box* (ROM + BiLSTM) y minimiza un costo de seguimiento, suavidad del actuador
y erosión del margen determinista (DSM).

**Qué falta.** La probabilidad de que el ACS falle no entra en la decisión. El margen
probabilístico (PSM) del árbol de eventos IS-MCDET se calcula *después* de elegir la acción.

**Qué se propone.** Primero volver tratable el problema: escribirlo de forma monolítica,
linealizar la dinámica del reactor en torno a puntos de operación y medir contra ese óptimo
la calidad de las metaheurísticas. Sobre esa base, formular un MPC estocástico (P1) cuyo
costo es un promedio ponderado sobre escenarios más una medida de riesgo. El resultado que
importa es **una decisión que cambia por causa del riesgo**, con intervalos de confianza y
en menos de 120 s por consulta. Las cotas certificadas del modelo no lineal pasan a ser una
línea secundaria.

---

## 2. Contexto mínimo

| Elemento | Qué es | Dónde está en el código |
|---|---|---|
| Planta | Modelo de alta resolución (HRM) del SMDFR que hace de reactor real | `SMDFR_hrm_MCDET_2022b_control_CMEX` |
| Accidente | LOCA de pequeña ruptura; el ACS se demanda sobre 1180 K; límite $U_y=1244$ K | `run_phase2c_full_loop.m` |
| Actuador | Bomba MHD: corriente $u$ [A] → caudal $\dot m_{\rm ACS}$ [kg/s] (modelo ECM de Nowak et al., 2022) | `CurrentToMassFlow.m` |
| Modelo de predicción | ROM físico + BiLSTM con MC-dropout que pronostica fuga $\hat d_k$ y error del ROM $\hat\Delta_k$ | `DET_2.m`, `python/mc_inference.py` |
| Incertidumbre discreta | ACS con 5 componentes; fallas y reparaciones generan ramas de un árbol dinámico de eventos | `DET_2.m` |
| Riesgo | PSM a partir del cuantil 95 % de la PCT máxima, ponderado por muestreo por importancia | `DET_scheduler_postprocessing_perc_online.m` |
| Optimizador | Luus–Jaakola (principal) y algoritmo genético; PID como referencia | `luus_jaakola_det.m`, `ga_det.m` |

Datos operativos: $\Delta t=0.5$ s, horizonte $N=25$ pasos (12.5 s), 61 consultas en
$t\in[1160,1190]$ s, 26 variables de decisión, cotas $u\in[56.69,\,56.90]$ A
(originalmente 56.75 A). La bomba arranca en 56.687 A y entrega ≈ 205 kg/s por A;
$u_{\min}$ todavía inyecta 0.58 kg/s y el ACS sólo puede enfriar.

---

## 3. Formulación del problema

### 3.1 Modelo de predicción

Para $k=0,\dots,N$ desde el instante de consulta $t_q$:

$$
\dot m_k = f_{\rm ACS}(u_k),\qquad
x_{k+1} = f(x_k,\sigma_k,\dot m_k,\hat d_k),\qquad
\hat y_k = h(x_k)+\hat\Delta_k,\qquad
\sigma_{k+1}\sim P(\cdot\mid\sigma_k,\hat y_k;\lambda,\mu).
$$

$x_k$ es el estado del ROM, $\hat y_k$ la PCT estimada y $\sigma_k$ el modo del ACS
(sin demanda / demandado y disponible / demandado y fallado). La incertidumbre tiene una
parte continua (pronóstico de fuga y error del ROM) y una discreta (instantes de falla y
reparación).

### 3.2 P0: el problema implementado

$$
\min_{u_0,\dots,u_N}\;
\underbrace{\sum_{j\in\mathcal I(u)}(\hat y_j-y_{\rm ref})^2}_{J_T}
+\underbrace{\lambda_m\sum_{k}\big[\,|\dot m_{k+1}-\dot m_k|-\Delta\dot m\,\big]_+^2}_{J_m}
+\underbrace{\lambda_s\sum_{j\in\mathcal I(u)}\big[\,\mathrm{DSM}_j-\tau\,\big]_+^2}_{J_s}
\quad\text{s.a.}\quad u_{\min}\le u_k\le u_{\max}
$$

con $\mathrm{DSM}_j=(\hat y_j-y_{\rm ref})/(U_y-y_{\rm ref})$, $y_{\rm ref}=1179.27$ K,
$\Delta\dot m=5$ kg/s, $\lambda_m=10^{-4}$, $\tau=0.00238$ y $\lambda_s\in\{0,10^4\}$
(`fitness_cost.m`). Se aplica sólo $u_0^\star$ y se repite en la consulta siguiente.

### 3.3 Las cuatro brechas de P0

| # | Brecha | Evidencia | Consecuencia |
|---|---|---|---|
| 1 | El costo no es una esperanza | `fitness_cost.m` concatena muestras de todos los nodos de `forest_data{1}` sin $p(\omega)$ | Candidatos con distinta topología de árbol se comparan sobre bases distintas |
| 2 | El riesgo llega tarde | El PSM se calcula en el posprocesamiento | El riesgo probabilístico no influye en la acción |
| 3 | Horizonte corto | En 12.5 s el árbol casi siempre queda sólo en la raíz | La falla del ACS es invisible para el optimizador |
| 4 | Margen de una trayectoria | $J_s$ usa el DSM de la trayectoria mediana | No mide probabilidad ni severidad de exceder $U_y$ |

Además: el chequeo `feasible` ($|\mathrm{DSM}|\le1$) no se impone en Luus–Jaakola pero sí
filtra en el algoritmo genético. Hay que unificar esa política antes de comparar optimizadores.

### 3.4 P1: el problema objetivo

$$
\begin{aligned}
\min_{\pi}\quad & \mathbb E_p\Big[\sum_{k=0}^{N}\ell(\hat y_k,\dot m_k)\Big]+\lambda_R\,\rho(Z^\pi),
\qquad Z^\pi=\max_{0\le k\le N_R}\hat y_k-U_y\\
\text{s.a.}\quad & \mathbb P_p\big(Z^\pi>0\mid\mathcal I_q\big)\le\varepsilon,\qquad
u_k=\pi_k(\mathcal I_k)\in[u_{\min},u_{\max}],\qquad\text{dinámica de §3.1.}
\end{aligned}
$$

- $\ell$ conserva los términos de seguimiento y suavidad de P0.
- $\mathbb E_p$ pondera cada ruta del árbol por su probabilidad (corrige la brecha 1).
- $\rho$ o la restricción con $\varepsilon$ meten el riesgo en cada candidato (brechas 2 y 4).
- $N_R\ge N$ separa el horizonte de riesgo del de seguimiento (brecha 3).
- $\pi$ es una política no anticipativa: sólo usa lo ya observado.

**Dos niveles.** *P1a*: una sola secuencia $u_{0:N}$ para todos los escenarios (misma
dimensión que P0; cambia sólo el costo). *P1b*: multietapa, con acciones distintas
después de observar éxito o falla del ACS y acciones idénticas antes del evento.
Se resuelve P1a primero; P1b sólo si P1a muestra que el riesgo cambia decisiones.

### 3.5 Medidas de riesgo candidatas

| Medida | Definición | Ventaja | Dificultad |
|---|---|---|---|
| Probabilidad de excedencia | $\mathbb P(Z>0)$ | Lenguaje regulatorio directo | Discontinua; eventos raros exigen IS |
| PSM por cuantil | $(U_y-q_{0.95}(\max\hat y))/(U_y-y_{\rm nom})$ | Continuidad con IS-MCDET | No convexo; ignora la cola más allá del cuantil |
| CVaR del exceso | $\mathbb E[Z\mid Z\ge\mathrm{VaR}_\alpha]$ | Mide severidad; convexo en $Z$ | Elegir $\alpha$; unidades en K |

El PSM y la probabilidad no son intercambiables. El evento, el horizonte de misión y
$\alpha$ o $\varepsilon$ se fijan **antes** de comparar métodos.

### 3.6 Aproximación computable

Con escenarios $\omega_i\sim q$, pesos $w_i=p(\omega_i)/q(\omega_i)$ calculados sobre la
trayectoria completa y $\bar w_i=w_i/\sum_l w_l$:

$$
\hat J(u)=\sum_i\bar w_i J(u,\omega_i),\qquad
\hat P(u)=\sum_i\bar w_i\,\mathbf 1\{Z_i(u)>0\},\qquad
\widehat{\mathrm{CVaR}}_\alpha(u)=\min_\eta\;\eta+\tfrac{1}{1-\alpha}\sum_i\bar w_i[Z_i(u)-\eta]_+ .
$$

Reglas: escenarios **comunes** para todos los candidatos; reportar
$\mathrm{ESS}=(\sum w_i)^2/\sum w_i^2$ e intervalos de confianza; elegir con unos
escenarios y evaluar con otros. Si las tasas de falla dependen de la acción (por ejemplo,
carga del componente), $p$ debe ser condicional al estado y a la acción.

---

## 4. Lo que ya sabemos

**Lazo cerrado** (61 consultas, Luus–Jaakola, semilla 1, mismo accidente;
`results/results_retuned_lambda1e4/`):

| Indicador | A: base ($\lambda_s=0$, $u_{\max}=56.75$) | D: reajustado ($\lambda_s=10^4$, $u_{\max}=56.90$) |
|---|---:|---:|
| PCT máxima / margen mínimo | 1181.21 K / 62.79 K | 1181.21 K / 62.79 K |
| IAE | 35.16 K·s | 28.65 K·s |
| Error en régimen | +0.07 K | +0.70 K |
| Variación total de $u$ | 0.86 A | 2.54 A |
| Tiempo mediano por consulta | 351.8 s | 76.0 s |

Se puede afirmar: la envolvente no cambia (el pico ocurre en la primera consulta), el
seguimiento mejora a costa de un actuador más activo y el cómputo es 4.6× más rápido. No
se puede afirmar que el controlador reduzca el riesgo, porque el riesgo no se evaluó dentro
de la decisión y hay una sola realización por configuración. D cambia cota y peso a la vez.

**Lección estructural** (`docs/DSM_PENALTY_RESULTS.md`): un término de seguridad sólo
desplaza la decisión cuando el óptimo de seguimiento es interior. En fase aguda seguimiento
y seguridad piden lo mismo (más refrigerante) y el óptimo satura en la cota: 0 de 49 ajustes
lo mueven. En convergencia, 2 de 49 lo mueven, y sólo con un umbral bajo la banda alcanzable.
Para observar el efecto del riesgo hacen falta escenarios donde seguimiento y seguridad
compitan.

**Base computacional.**

- *Reproducibilidad*: el piloto en R2024b da $J=0.362026$ frente a $0.326151$ archivado
  en R2024a (+11.0 %). Causa no aislada ([`pilot_reproducibility.md`](pilot_reproducibility.md)).
- *Búsqueda*: Luus–Jaakola devolvía acciones hasta ~190× peores que una malla; se corrigió
  (techo geométrico del radio y parada relativa).
- *Tiempo de cálculo*: mediana de 76 s por consulta, dentro de la ventana de 120 s que la
  tesis adopta para un monitor de riesgo en línea. Un PSM por candidato con el criterio de
  Wilks (95/95, 59 corridas) multiplicaría ese costo por ≈ 59.
- *Cotas*: `m1_branch_and_bound.m` cierra la brecha en un ejemplo híbrido
  ($0.499024\le J^\star\le0.5$), pero para el gemelo real la única cota inferior es $J\ge0$
  ([`m1_cost_bound_results.md`](m1_cost_bound_results.md)).

---

## 5. Vías de exploración

Criterios de prioridad: ¿responde la pregunta de investigación? ¿Es factible con el código
actual? ¿Habilita otras vías?

### Lo que cambió la reunión del 22-10-2026

Diagnóstico compartido: ni siquiera el MPC determinista se resuelve con holgura, y meter
el riesgo probabilístico lo multiplica. **Antes de volverlo estocástico hay que volverlo
tratable.** De ahí salen tres vías nuevas que pasan al frente (V1–V3, propuestas por
R. Mena y H. Rozas) y una cuarta que se explora después (V6, propuesta por H. Rozas).
Acuerdo: buscar modelos linealizados del reactor (funciones de transferencia en torno a
puntos de operación), reformular el plan y volver a reunirse en 3–4 semanas.

### V1. Formulación monolítica — prioridad alta

- **Idea**: escribir P0 como un problema de optimización estándar en vez de hacerlo en la
  forma que exige la metaheurística (todo dentro del costo):

  $$
  \begin{aligned}
  \min_{\dot m,x,s,r}\ & \textstyle\sum_k(\hat y_k-y_{\rm ref})^2+\lambda_m\sum_k s_k^2+\lambda_s\sum_k r_k^2\\
  \text{s.a.}\ & x_{k+1}=f(x_k,\dot m_k,\hat d_k),\quad \hat y_k=h(x_k)+\hat\Delta_k\\
  & |\dot m_{k+1}-\dot m_k|\le\Delta\dot m+s_k,\quad \mathrm{DSM}_k\le\tau+r_k,\quad s_k,r_k\ge0\\
  & \dot m_{\min}\le\dot m_k\le\dot m_{\max}.
  \end{aligned}
  $$

  Las bisagras $[z]_+^2$ de P0 se reformulan **exactamente** con holguras. El margen puede
  pasar a una restricción dura y lineal: $\hat y_k\le y_{\rm ref}+\tau(U_y-y_{\rm ref})$.
  Cuadrados y valores absolutos son manejables por solvers actuales, con aproximaciones
  cónicas si hace falta.
- **Hallazgo que la respalda**: en $u\in[56.69,\,56.90]$ A el mapa de la bomba es afín,
  $\dot m\approx204.6\,(u-56.687)$, con error menor que 0.007 kg/s (0.015 % del rango).
  Se puede decidir el caudal directamente y recuperar $u$. Lo único no lineal que queda es
  la dinámica térmica del ROM ($f$, $h$), con el retardo de la transferencia de calor.

### V2. Modelo lineal del reactor — prioridad alta

- **Idea**: identificar
  $\delta x_{k+1}=A\,\delta x_k+B\,\delta\dot m_k+E\,\delta\hat d_k$,
  $\hat y_k=\bar y+C\,\delta x_k+\hat\Delta_k$ con respuestas del ROM (escalones o señales
  pseudoaleatorias) en torno al estado de cada consulta, o en varios puntos de operación.
  El retardo térmico entra como estados adicionales. Con V1, el problema completo es un
  **QP convexo** (MPC lineal-cuadrático clásico).
- **Qué abre**: MPC robusto (tubos), distribucionalmente robusto y contextual; restricción
  probabilística con error gaussiano como ajuste determinista,
  $\mathbb P(\hat y_k>U_y)\le\varepsilon\iff\bar y_k+z_{1-\varepsilon}\sigma_k\le U_y$, con
  $\sigma_k$ estimada de las 59 trayectorias de MC-dropout; y un modelo por modo del ACS
  (sistema conmutado) para las fallas.
- **Criterio de validez**: la autoridad del control sobre la PCT en el horizonte es ≈ 1 K
  (0.975 K en $t=1180$ s con la cota original), así que el error de predicción lineal debe
  ser mucho menor que 1 K. En la fase aguda puede no bastar un solo punto de operación.
- **Primer paso**: buscar en la literatura linealizaciones o funciones de transferencia del
  reactor; si no hay, identificarlas desde el ROM.

### V3. Referencia exacta para las metaheurísticas — prioridad alta

- **Idea**: resolver los mismos estados de consulta con Luus–Jaakola, el algoritmo genético
  y el solver de V1–V2; re-simular las acciones en el ROM; medir brecha, acción y tiempo.
- **Dos brechas**: sobre el modelo lineal (exacta, con la cota dual del solver) y sobre el
  ROM (la que importa físicamente).
- **Por qué**: da la referencia comparativa que las cotas del modelo no lineal (V10) no
  lograron dar; hoy sólo se sabe que $J\ge0$.

### V4. Riesgo en cada candidato — prioridad alta

- **Idea**: sobre la base V1–V2, primero la restricción probabilística ajustada por cuantil;
  luego P1a con $\hat J+\lambda_R\hat\rho$ (o $\hat P\le\varepsilon$) sobre escenarios comunes
  con pesos IS.
- **Primer experimento**: en un instante fijo, $M$ escenarios (pronóstico × eventos del ACS);
  calcular $\hat J$, $\hat P$ y $\widehat{\mathrm{CVaR}}$ con los mismos escenarios para todo
  candidato; comparar la acción de P0 y P1a para varios $\lambda_R$ o $\varepsilon$ y
  re-simularlas en el ROM.
- **Riesgo**: que en los casos actuales P1a coincida con P0 (por la lección estructural).
  Por eso V4 va con V5 y con casos donde seguimiento y seguridad compitan.

### V5. Horizonte de riesgo extendido — prioridad alta

- **Idea**: evaluar $\rho$ sobre $N_R\gg N$ (hasta el tiempo de misión), con una política
  terminal simple más allá de $N$ (por ejemplo, mantener $u_N$ o la ley de procedimiento).
- **Por qué**: con 12.5 s el árbol casi no ramifica y la falla del ACS no se ve.

### V6. Escenario orientado a la decisión — prioridad media

- **Idea**: fuera de línea, resolver el MPC estocástico con muchos escenarios para muchos
  contextos $c$ (estado, pronóstico) y entrenar $g_\theta:c\mapsto\hat\xi$ tal que el MPC
  determinista alimentado con el escenario único $\hat\xi$ dé una decisión con garantías
  similares a las del problema estocástico. El predictor no minimiza el error de pronóstico,
  sino la calidad de la decisión. En línea sólo se resuelve un MPC determinista.
- **Condiciones**: es una idea conceptual; no se sabe si funciona aquí y las garantías hay
  que demostrarlas. Exige resolver miles de MPC fuera de línea, lo que requiere V1–V2.
- **Base teórica**: Homem-de-Mello, Valencia, F. Lagos y G. Lagos (2026), *Forecasting
  outside the Box*, Math. Program., doi:10.1007/s10107-026-02359-9 (enviado por H. Rozas).
  Su teorema exige dos etapas, segunda etapa lineal y recurso fijo; el análisis de
  aplicabilidad está en [`monolithic_formulation_plan.md`](monolithic_formulation_plan.md) §7.
- **Desarrollo de V1–V3**: la escalera de modelos N0–N3, la formulación monolítica y el
  plan F0–F5 están en [`monolithic_formulation_plan.md`](monolithic_formulation_plan.md).

### V7. Cómputo del modelo no lineal — prioridad media

- **Idea**: bloqueo de movimientos $u_k=\sum_{j=1}^{n_b}\theta_j\phi_j(k)$ con
  $n_b\approx3$–5 tramos constantes; arranque en caliente; presupuesto fijo de evaluaciones.
- **Cuándo hace falta**: si el modelo lineal no es válido en el transitorio y el riesgo debe
  evaluarse sobre el ROM. P1a multiplica cada evaluación por $M$ ($M=59$ para un PSM de
  Wilks); sin V7 no cabe en 120 s.
- **Hipótesis a verificar**: la respuesta térmica es lenta frente a 0.5 s y se pierde poco
  desempeño.

### V8. MPC multietapa (P1b) — prioridad media

- **Idea**: árbol de decisión con recurso: acción común hasta el evento del ACS y acciones
  $u^{\rm ok}$, $u^{\rm falla}$ después, con no anticipación explícita.
- **Condición de entrada**: que V4 haya mostrado al menos una decisión sensible al riesgo.

### V9. Varias protecciones y sistema pasivo — prioridad media

- **Idea**: extender $\sigma_k$ a dos sistemas activos (disponible, demandado, activado,
  fallado, reparado; retardos; recursos compartidos) y al tapón de sal pasivo, con demanda y
  disponibilidad propias y dependencia de la secuencia activa.
- **Por qué importa**: con un solo sistema activo, una falla del ACS no admite respaldo y el
  riesgo sólo puede actuar de forma preventiva. Con dos aparece la decisión de qué activar
  primero y qué hacer si uno falla; A→B y B→A pueden producir riesgos distintos.
- **Requisito**: poda del árbol con regla auditable y cota del riesgo omitido.

### V10. Cotas del modelo no lineal — prioridad baja, línea secundaria

- **Idea**: encierros del ROM/DET por cajas de control (aritmética intervalar o relajación
  validada) para cotas inferiores reales y una brecha $U-L$ con significado.
- **Estado**: algoritmo y ejemplo híbrido listos; falta el encierro térmico del ROM. Si no
  es viable, V3 da la referencia práctica.
- **Actualización (9-11-2026)**: con la réplica algebraica en Julia (F1), cada evaluación
  pasa de minutos a milisegundos, y la vía se vuelve viable para el ROM real. Queda
  incorporada como etapa **F3c** del plan monolítico
  ([`monolithic_formulation_plan.md`](monolithic_formulation_plan.md) §8.5): certificado
  local, NLP sobre `ode1be`, evidencia global y cota inferior rigurosa sobre la caja de
  control.

---

## 6. Plan por etapas

```
E0 Línea base ─► E1 Monolítico ─► E2 Referencia ─► E3 Riesgo en la ─► E4 Lazo cerrado ─► E5 Multietapa y
   congelada        y lineal         exacta           decisión           y tiempo real      protecciones
                    (V1–V2)          (V3)             (V4–V5)            (V6 / V7)          (V8–V9)
                 └─ V10 (cotas) secundaria ─┘
```

| Etapa | Entregable | Criterio de cierre |
|---|---|---|
| **E0** Línea base congelada | Archivo de configuración único (modelo, estado inicial, horizonte, semillas, factibilidad); el mismo vector $u$ evaluado en R2024a y R2024b | Diferencia del 11 % explicada o acotada; un comando reproducible |
| **E1** Monolítico y lineal | Modelo lineal identificado desde el ROM en al menos dos puntos (fase aguda y convergencia); P0 escrito como QP con holguras | Error de predicción lineal mucho menor que 1 K en el horizonte; el QP reproduce P0 en puntos de prueba |
| **E2** Referencia exacta | Luus–Jaakola, genético y solver sobre los mismos estados; acciones re-simuladas en el ROM | Brecha de optimalidad y tiempo de las metaheurísticas medidos frente a la referencia |
| **E3** Riesgo en la decisión | Restricción probabilística y escenarios sobre la base lineal; envoltorio de `DET_2` que devuelve una tabla por escenario (ruta, $p$, $q$, $\hat y_k$, modo) | Al menos un caso con $u^\star_{\rm P1a}\ne u^\star_{\rm P0}$, diferencia de riesgo fuera de los intervalos de confianza y costo de seguimiento declarado |
| **E4** Lazo cerrado y tiempo real | Varias realizaciones con semillas independientes; P0, P1a y PID con igual presupuesto; V6 como alternativa en línea | Frontera costo–riesgo fuera de muestra; menos de 120 s por consulta; un caso de falla documentado |
| **E5** Extensiones | Modos para dos sistemas activos y pasivo; política P1b | Casos sintéticos A→B, B→A, falla común, ACS fallido y activación pasiva, con trazas causales |

**Casos de prueba para E2–E4.** $t=1161$ s (fase aguda, control saturado); $t=1180$ s
(convergencia, óptimo interior); fuga severa $B=15$ (`generate_severe_scenario.m`); ACS
degradado (`vel_ACS` reducido en el HRM con el mismo script; la serie archivada
`LJ_test_sampled_vel_ACS_09.mat`, con $t_q=1300$ s, no se distribuye aquí); y una falla
del ACS forzada dentro del horizonte.

---

## 7. Decisiones abiertas

1. **Puntos de operación**: uno o varios modelos lineales; cuánto error lineal se tolera en la
   fase aguda.
2. **Solver y forma del problema**: QP, cónico o no lineal; comercial o abierto; margen como
   restricción dura o suave.
3. **Medida y evento**: probabilidad, PSM por cuantil o CVaR; valor de $\varepsilon$ o $\alpha$;
   horizonte de misión.
4. **Costo o restricción**: $\lambda_R\rho$ es fácil de optimizar pero exige calibrar
   $\lambda_R$; $\mathbb P\le\varepsilon$ es interpretable pero puede volver infactible el problema.
5. **Presupuesto**: 120 s por consulta; cuántos escenarios caben y si hace falta el escenario
   aprendido (V6).
6. **Segunda protección**: qué sistema activo adicional, con qué tasas, retardos y dependencias
   con el tapón pasivo.

---

## 8. Comandos de verificación

```text
python python/verify_fit_poly_coefficients.py
matlab -batch "addpath('matlab'); test_certified_cost_lower_bound_box; test_m1_branch_and_bound"
matlab -batch "addpath('matlab'); m1_pilot_audit"
matlab -batch "addpath('matlab'); addpath('models'); pilot_single_run"
```

El último requiere MATLAB R2024a, Simulink, modelos y MEX; escribe `pilot_result.mat` en el
directorio de ejecución, así que conviene correrlo desde un directorio de resultados dedicado.

---

## 9. Referencias metodológicas y límites de transferencia

- Rockafellar y Uryasev (2000), *Optimization of Conditional Value-at-Risk*,
  https://doi.org/10.21314/jor.2000.038. Base de la representación optimizable de CVaR; se
  propone como alternativa de cola, no como sustituto automático del PSM.
- Shapiro (2003), *Monte Carlo Sampling Methods*, https://doi.org/10.1016/S0927-0507(03)10006-0.
  SAA y validación estadística del valor objetivo; sus estimadores no reemplazan una cota
  determinista.
- Campi y Garatti (2008), *The Exact Feasibility of Randomized Solutions of Uncertain Convex
  Programs*, https://doi.org/10.1137/07069821X. Garantías por escenarios bajo convexidad; el
  gemelo híbrido con modos discretos no las satisface automáticamente.
