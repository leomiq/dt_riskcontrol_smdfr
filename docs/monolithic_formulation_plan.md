# Hacia una formulación monolítica del MPC del SMDFR

Versión: 22 de octubre de 2026. Documento compañero de
[`presentations/monolithic_formulation.pdf`](../presentations/monolithic_formulation.pdf).
Desarrolla las vías V1 (formulación monolítica), V2 (modelo lineal) y V3 (referencia exacta)
de [`implementation_route.md`](implementation_route.md), y conecta la vía V6 con el artículo
que propuso H. Rozas.

---

## 1. Para qué sirve este documento

La reunión del 22-10-2026 concluyó que, antes de volver estocástico el MPC, hay que volverlo
tratable. Después aparecieron cuatro caminos posibles:

1. linealizar por pequeña perturbación y obtener funciones de transferencia, como en el
   HPR1000;
2. usar una respuesta estacionaria, como la de `SA_DualFluid.m`;
3. escribir las ecuaciones de cinética puntual (PKE) y termohidráulica (TH) del SMDFR como
   restricciones, fusionando el simulador con el problema de control;
4. reutilizar los modelos en Julia del secundamiento en ARAMIS.

Este documento ordena esos caminos. La tesis central es que **no son alternativas, sino
niveles de un mismo modelo**. Si las ecuaciones del ROM se escriben una sola vez en forma
algebraica, de ellas salen el modelo estacionario, el linealizado en un punto, el linealizado
a lo largo de la trayectoria y el modelo no lineal completo. El problema monolítico se
construye sobre esa única fuente.

---

## 2. Inventario del material

| Material | Qué aporta | Cómo se usa aquí | Límite |
|---|---|---|---|
| Homem-de-Mello, Valencia, F. Lagos y G. Lagos (2026), *Forecasting outside the Box*, Math. Program., doi:10.1007/s10107-026-02359-9 (propuesto por H. Rozas) | En programas estocásticos de dos etapas con recurso fijo, costos fijos y segunda etapa lineal (FRFC), existe un **escenario óptimo único**: resolver el problema con ese escenario da la solución estocástica. Con información contextual, se aprende un pronóstico puntual por optimización binivel. | Base teórica de V6 (escenario orientado a la decisión), §7 | Requiere segunda etapa **lineal** con recurso fijo; nuestro MPC es no lineal, multietapa y con costo cuadrático |
| Zhang, Sun, Qiu, Pu y Wei (2022), *Transfer function modeling and simulation of HPR1000*, Ann. Nucl. Energy 166:108689 | Procedimiento: modelo no lineal por conservación → estado estacionario → linealización por pequeña perturbación → espacio de estados $\dot x=Ax+Bu$ → $G(s)=C(sI-A)^{-1}B+D$. Validado contra el no lineal con perturbaciones típicas | Nivel N1 (§4) | Un punto de operación; las perturbaciones de validación son pequeñas (50 pcm, 5 % de válvula) |
| Wang, He, Macián-Juan y Wang (2021), *One-dimensional transient analysis of the dual-fluid reactor system*, Ann. Nucl. Energy 162:108481 | PKE con **precursores circulantes** (pérdida en el núcleo y reingreso tras el lazo) acoplada a balances de calor nodalizados; datos del DFR-250M (Λ, β_i, λ_i, α_f, α_l, tiempos de tránsito) | Referencia de las ecuaciones de N3 y de su validación física | Analiza transitorios abiertos, no control |
| `SA_DualFluid.m` y *Presentazione AI* (G. Masotti y R. Mauri) | Modelo TH **estacionario** del núcleo en 12 nodos: un sistema lineal 5×5 por nodo. Su resultado (`models/SA_DualFluid.mat`) **inicializa el HRM y el ROM**. Su análisis de sensibilidad muestra que $h_{fw}$ domina la incertidumbre de $T_f$ (η² ≈ 0.86) | Nivel N0 (§4) y estado inicial de N1–N3 | Sin dinámica ni cinética |
| HRM y ROM del repositorio (`models/*.slx`) | Ecuaciones PKE + TH ya implementadas: 12 nodos (HRM) y 0-D (ROM) | **Fuente de las restricciones** (§3) | Están en bloques Simulink; hay que extraerlas y validarlas |
| Modelos Julia del secundamiento ARAMIS (`leomiq/smdfr_control`) | Mismo reactor y parámetros (PKE de 6 grupos + TH de 3 nodos, como DAE) y **un MPC monolítico no lineal ya escrito** en JuMP + Ipopt con Euler implícito | Plantilla de N3 (§6) | El control es la reactividad externa, no el caudal del ACS; no hay salida de PCT ni escalado |

---

## 3. El modelo de predicción, escrito completo

Las ecuaciones se extrajeron del archivo `models/rom_MCDET_v2_2022b_control_CMEX.slx`
(bloques MATLAB Function y cableado) con los parámetros de
`matlab/parameters_dualfluid_rom_v2.m`. El ROM es **de dimensión cero**: un nodo de
combustible, uno de pared y uno de refrigerante, más la cinética.

**Cinética puntual con precursores circulantes** (6 grupos):

$$
\dot N=\frac{\rho-\beta}{\Lambda}N+\sum_i\lambda_iC_i,\qquad
\dot C_i=\frac{\beta_i}{\Lambda}N-(\lambda_i+\lambda_c)C_i+\lambda_eC^e_i,\qquad
\dot C^e_i=\lambda_cC_i-(\lambda_i+\lambda_e)C^e_i,
$$

$$
\rho=\rho_0+\alpha_f(T_f-\bar T_{f0})+\alpha_c(T_c-\bar T_{c0}),\qquad P=\tfrac{P_0}{N_0}N .
$$

$\lambda_c=1/\tau_c$ y $\lambda_e=1/\tau_e$ representan la salida de precursores del núcleo y
su reingreso desde el lazo externo. **Comprobación**: con estos parámetros, la reactividad
estacionaria que exige la pérdida de precursores es
$\beta-\sum_i\beta_i\lambda_i/(\lambda_i+\lambda_c-\lambda_c\lambda_e/(\lambda_i+\lambda_e))
=1.1105\cdot10^{-3}$, que coincide con $\rho_0=1.11\cdot10^{-3}$ del modelo.

**Termohidráulica** (nodo promedio, $T_e=2T-T_{in}$):

$$
M_fc_{pf}\dot T_f=P+\dot m_fc_{pf}(T_{in,f}-T_{e,f})-h_{fw}A_{fw}(T_f-T_w)
$$

$$
M_wc_{pw}\dot T_w=h_{fw}A_{fw}(T_f-T_w)-h_{wc}(\dot m_c)A_{wc}(T_w-T_c)
$$

$$
M_cc_{pc}\dot T_c=\dot m_cc_{pc}(T_{in,c}-T_{e,c})+h_{wc}(\dot m_c)A_{wc}(T_w-T_c)
$$

**Inventario de refrigerante y actuador**: el control actúa a través del inventario.

$$
\dot M=\dot m_{\rm leak}(t)+\sigma\,\dot m_{\rm ACS}(u),\qquad
\dot m_c=\dot m_{c0}\frac{M}{M_0},\qquad
h_{wc}(\dot m_c)=h_{wc0}\Big(\frac{\dot m_c}{5550}\Big)^{0.8},
$$

con $\dot m_{\rm ACS}(u)\approx204.6\,(u-56.687)$ en la ventana de operación y $\sigma\in\{0,1\}$
el estado del ACS. En el HRM, la fuga sigue $\dot m_{\rm leak}=-k_{\rm leak}+A(t-t_f)^{1/B}$;
en el MPC se usa el pronóstico de la BiLSTM.

**Salida (PCT)**:

$$
\hat y=T_w+p_4(\dot m_{\rm net}+c)+0.4,\qquad \dot m_{\rm net}=\dot m_{\rm leak}+\dot m_{\rm ACS},
$$

donde $p_4$ es el polinomio de grado 4 de `fit_poly` (diferencia HRM−ROM ajustada contra el
caudal neto).

**Dimensión**: 16 estados ($N$, 6 $C_i$, 6 $C^e_i$, $T_f$, $T_w$, $T_c$, $M$) y una
entrada de control.

### 3.1 Puntos a confirmar antes de escribir restricciones

1. **La corrección BiLSTM de la PCT no entra al ROM.** La entrada `pred_delta_tw_max` llega
   al subsistema `SMDFR 1`, pero dentro no está conectada. La corrección efectiva es $p_4$.
   La BiLSTM sólo aporta el pronóstico de la fuga. Hay que confirmar si es intencional y
   corregir la descripción del gemelo (README, presentación anterior y tesis hablan de
   $\hat\Delta_k$ como salida de la red).
2. `Tbar_f0 = mean(Te_f0, Tin_f0)` usa el segundo argumento como dimensión, así que
   $\bar T_{f0}=T_{e,f0}$ (y lo mismo para el refrigerante). Es inocuo si es lo pretendido.
3. La temperatura de entrada del refrigerante usa una relación empírica
   $T_{in,c}=T_{e,c0}-Q(\dot m_c)/(\dot m_cc_{pc})$ con
   $Q=Q_0\dot m_c/2775-(5500-\dot m_c)$ y $Q_0=0.578\cdot8.2564\cdot10^7$ W. Conviene
   documentar su origen.
4. Hay un desplazamiento de −2.5 K en el término constante de $p_4$
   ([`pilot_reproducibility.md`](pilot_reproducibility.md)) además del +0.4 K de la salida.

---

## 4. Escalera de modelos: un modelo, cuatro usos

| Nivel | Modelo | Cómo se obtiene | Problema de control | Qué captura | Qué pierde |
|---|---|---|---|---|---|
| **N0** | Estacionario | Resolver §3 con derivadas nulas (como `SA_DualFluid.m`) | Mapa estático $\dot m_c\mapsto T_w$; restricción terminal o chequeo de autoridad | Ganancia en régimen | Toda la dinámica: inventario, precursores, retardo |
| **N1** | Lineal invariante (LTI) | Pequeña perturbación en el estado de la consulta: $A=\partial f/\partial x$, $B=\partial f/\partial u$, $E=\partial f/\partial d$ (método HPR1000); $G(s)=C(sI-A)^{-1}B+D$ | QP convexo; MPC lineal clásico, robusto o con restricción probabilística | Dinámica cerca del punto | Validez lejos del punto (fase aguda) |
| **N2** | Lineal variante (LTV) | Linealizar a lo largo de la trayectoria predicha en la consulta anterior (iteración en tiempo real, SQP) | Un QP por iteración | Transitorio grande, si se itera | Garantía de convergencia en una sola iteración |
| **N3** | No lineal completo | Ecuaciones de §3 discretizadas como restricciones (colocación o disparo múltiple) | NLP (Ipopt) | Todo el ROM | Convexidad; tiempo de cálculo mayor |

**Recomendación**: escribir el ROM **una sola vez** en un lenguaje algebraico con
diferenciación automática (JuMP, CasADi o ModelingToolkit). N0 es su solución estacionaria;
N1 y N2 son sus jacobianos, en un punto o a lo largo de una trayectoria; N3 es el modelo
mismo. Así los cuatro niveles son consistentes por construcción, y N3 sirve de referencia
exacta para medir a N1, N2 y a las metaheurísticas (V3).

### 4.1 Escalas de tiempo y discretización

Las magnitudes salen de los parámetros del ROM:

| Fenómeno | Escala | Consecuencia con $\Delta t=0.5$ s |
|---|---|---|
| Neutrones prontos | $\Lambda/(\beta-\rho_0)\approx5.5\cdot10^{-4}$ s | Euler explícito es inestable; usar **aproximación de salto pronto** (PJA: $N$ algebraico) o integración implícita |
| Pared | $\tau_w\approx0.07$ s | Casi estática frente a $\Delta t$: candidata a ecuación algebraica |
| Refrigerante, combustible | $\tau_c\approx0.25$ s, $\tau_f\approx0.72$ s | Del orden de $\Delta t$: Euler implícito o subpasos |
| Precursores | 0.14–79 s; tránsitos de 5.4 s (núcleo) y 10 s (lazo) | Lentos; bien resueltos |
| Inventario | 10 kg/s durante 12.5 s cambian $M$ en 1.5 % y $\dot m_c$ en ≈ 85 kg/s | Integrador: el efecto del ACS se **acumula** en el horizonte |

La razón de rigidez entre neutrones prontos y el precursor más lento es del orden de
$10^5$. Con PJA y pared cuasiestática, el sistema dinámico efectivo tiene 14 estados.

### 4.2 Las no linealidades y su tratamiento

| Término | Tipo | En N1/N2 | En N3 |
|---|---|---|---|
| $\rho N$ | Bilineal | Jacobiano | Exacto |
| $\dot m_c(T_{in}-T_e)$ | Bilineal (control × estado) | Jacobiano | Exacto; o McCormick si se busca una relajación convexa |
| $h_{wc}\propto\dot m_c^{0.8}$ | Potencia | Jacobiano | Exacto |
| $Q(\dot m_c)/\dot m_c$ | Racional | Jacobiano | Exacto ($\dot m_c$ acotado lejos de 0) |
| $p_4(\dot m_{\rm net})$ | Polinomio de grado 4 | Jacobiano | Exacto |
| $f_{\rm ACS}(u)$ | Afín en la ventana | Exacto | Exacto; decidir $\dot m_{\rm ACS}$ directamente |
| $[z]_+^2$, $\lvert\Delta\dot m\rvert$ | Bisagras | Holguras (exacto) | Holguras (exacto) |

---

## 5. La formulación monolítica

Discretización en $K$ pasos de $\Delta t=0.5$ s (o subpasos), con Euler implícito o
colocación de Radau:

$$
\begin{aligned}
\min_{\dot m_{\rm ACS},\,x,\,s,\,r}\quad & \sum_{k=1}^{K}(\hat y_k-y_{\rm ref})^2+\lambda_m\sum_k s_k^2+\lambda_s\sum_k r_k^2\\
\text{s.a.}\quad & x_{k+1}=x_k+\Delta t\,f(x_{k+1},\dot m_{{\rm ACS},k},\hat d_k)\quad\text{(ROM de §3, implícito)}\\
& 0=g(x_k)\quad\text{(PJA y pared cuasiestática, opcional)}\\
& \hat y_k=T_{w,k}+p_4(\hat d_k+\dot m_{{\rm ACS},k}+c)+0.4\\
& \lvert\dot m_{{\rm ACS},k+1}-\dot m_{{\rm ACS},k}\rvert\le\Delta\dot m+s_k,\quad s_k\ge0\\
& \mathrm{DSM}_k\le\tau+r_k,\quad r_k\ge0\qquad(\text{o la versión dura }\hat y_k\le y_{\rm ref}+\tau(U_y-y_{\rm ref}))\\
& \dot m_{\min}\le\dot m_{{\rm ACS},k}\le\dot m_{\max},\qquad x_0=\hat x(t_q).
\end{aligned}
$$

**Tamaño**: con $K=25$, son 16 × 26 ≈ 420 variables de estado, 26 de control y unas 50
holguras. Es un NLP pequeño. Que se resuelva en segundos es una **hipótesis a verificar**;
la referencia actual es una mediana de 76 s con Luus–Jaakola.

**Ventajas frente a P0**: el costo y las restricciones se evalúan sobre un modelo explícito
(sin una caja negra de Simulink), el solver entrega multiplicadores y, en N1, una cota dual.
Además, el mismo problema admite capas de riesgo (§7).

---

## 6. Modelos Julia del secundamiento ARAMIS

Ubicación: `Secondments/Industrial secondment - ARAMIS/July 2023 - Research Progress/smdfr_control/smdfr_control/toy example/`
(repositorio `leomiq/smdfr_control`). Hay 28 archivos `.jl`, sin `Project.toml`; el código
usa sintaxis antigua de JuMP (`@NLconstraint`) y Julia 1.8.

**Qué hay**

- **El mismo reactor y los mismos parámetros**: PKE de 6 grupos con precursores circulantes y
  TH de 3 nodos (pleno inferior, núcleo y pleno superior), portados línea a línea de
  `parameters_dualfluid_rom.m` (`parameters_dualfluid_rom.jl`). Hay un juego de 12 nodos
  comentado.
- **Ecuaciones en forma DAE semiexplícita con matriz de masa** (`reactor_ode+rod_insertion.jl`,
  `function_odes!`): 23 estados diferenciales y 8 algebraicos (reactividad, potencia y
  temperaturas de salida).
- **Un MPC monolítico no lineal ya escrito**: `smdfr_mpc_test_4.jl`,
  `optimizer_implicit_euler` (líneas 200–253). Discretiza la DAE con Euler implícito como
  restricciones de JuMP, impone el término bilineal $\rho N/\Lambda$ con `@NLconstraint`
  (línea 220), lo resuelve con Ipopt y lo usa en lazo cerrado dentro de un `PeriodicCallback`.
  **Es la plantilla natural de N3.**
- **Linealización**: jacobianos con ForwardDiff (`smdfr_mpc_test_3.jl`) y una PKE linealizada
  a mano (`smdfr_linearized_odes.jl`); no hay funciones de transferencia.

**Qué falta o hay que corregir**

1. **El control es otro**: la entrada es la tasa de reactividad externa (±15 pcm/s) y
   $\dot m_c$ es un parámetro fijo. Hay que reemplazarlo por el inventario $M$ y el caudal del
   ACS de §3. Eso introduce el bilineal $\dot m_c(T_{in}-T_e)$ y la dependencia $h_{wc}(\dot m_c)$.
2. **La salida es otra**: se sigue la potencia, no la PCT. Hay que añadir
   $\hat y=T_w+p_4(\cdot)$ (o el máximo nodal con restricciones de epígrafe si se usan 12 nodos).
3. **Escalado**: el código trabaja en unidades absolutas ($N\sim10^{13}$, $C_i\sim10^{16}$,
   $P\sim10^8$ W). Para Ipopt hay que normalizar ($n=N/N_0$, $c_i=C_i/C_{i0}$, temperaturas
   como desviaciones).
4. **Estado inicial consistente**: $\bar T_{f0}=1288$ y $\bar T_{c0}=1060$ K se ajustaron por
   inspección. Hay que resolver $f(x,u)=0$ (N0) y partir de ese equilibrio, o del estado
   del HRM en $t_q$.
5. **Linealización de `smdfr_mpc_test_3.jl`**: ignora la matriz de masa y omite el término
   afín $f(x_0)-Ax_0$. No debe reutilizarse tal cual.
6. **Higiene**: variables globales, cadenas de `include` y extracción manual de coeficientes
   (`Symbolics.coeff` con parches a mano). Conviene escribir el modelo directamente en JuMP
   moderno, con los parámetros en una estructura y un entorno fijado.

**Valoración**: la física necesaria ya está escrita en Julia, y el esqueleto de N3 existe.
El trabajo es de **adaptación**: cambiar la entrada al caudal del ACS, la salida a la PCT,
escalar y fijar el estado inicial. Hay dos caminos:

- **Julia** (JuMP + Ipopt): reutiliza la plantilla y es rápido; habría que reproducir la
  BiLSTM o pasar sus pronósticos como datos.
- **CasADi desde MATLAB o Python**: queda más cerca del lazo actual (Simulink y BiLSTM),
  pero hay que reescribir el modelo.

La recomendación es Julia para F1–F3 (prototipo y referencia exacta, fuera del lazo) y
decidir el lenguaje de integración en F5.

---

## 7. El riesgo sobre la base monolítica

Hay tres capas, de menor a mayor costo:

1. **Restricción probabilística analítica (N1/N2)**: si el error de pronóstico es gaussiano,
   $\mathbb P(\hat y_k>U_y)\le\varepsilon$ equivale a $\bar y_k+z_{1-\varepsilon}\sigma_k\le U_y$,
   con $\sigma_k$ obtenida de las 59 trayectorias de MC-dropout propagadas por el modelo
   lineal. Sigue siendo un QP.
2. **Escenarios (SAA)**: $M$ trayectorias de fuga y eventos del ACS. Con N1/N2 es un QP de
   tamaño $M$ veces mayor; con N3, un NLP con estructura de bloques.
3. **Escenario óptimo** (Homem-de-Mello et al., 2026): resolver un **único** escenario que
   reproduzca la decisión estocástica.

### 7.1 ¿Aplica el teorema del escenario óptimo?

| Supuesto del artículo | En nuestro MPC | Cómo cumplirlo |
|---|---|---|
| Dos etapas: decisión ahora, recurso después | MPC multietapa | Aproximación estándar: $u_0$ como primera etapa y el resto del horizonte como recurso |
| Segunda etapa **lineal** (LP) | Costo cuadrático, dinámica no lineal | Modelo N1/N2 y costo de seguimiento en norma $L_1$ ($y^+,y^-$ como en el *newsvendor*) |
| **Recurso fijo** ($W$, $q$ deterministas) | Si la incertidumbre es la fuga o el error del modelo, entra en el lado derecho $h$ ✓ | Cumple para el pronóstico de fuga |
| | La falla del ACS multiplica el efecto de controles futuros, así que vuelve aleatoria $W$ ✗ | Sólo admisible si la falla afecta a la primera etapa ($T$ aleatoria, caso ii del teorema); si no, queda fuera del teorema |
| Recurso relativamente completo (supuesto 2: duales acotados) | Las restricciones suaves con holguras penalizadas garantizan $Q$ finita | ✓ por diseño |

**Lectura**: el teorema cubre la incertidumbre de **pronóstico** (fuga, error del ROM) con un
MPC lineal y costo $L_1$. No cubre directamente la falla del ACS. Para ella sirven las
capas 1–2 y el MPC multietapa (V8).

### 7.2 Experimento mínimo inspirado en el artículo

En el *newsvendor*, el escenario óptimo es un **cuantil** de la demanda. El código ya
calcula los percentiles 50 y 95 del ensamble de MC-dropout (`Y_denorm_50`, `Y_denorm_95`),
y P0 usa la mediana. El experimento es: para varios niveles $q$, resolver el MPC determinista
con el percentil $q$ de la fuga y comparar su acción y su costo con los de la solución SAA,
que usa las 59 trayectorias. Si algún $q$ reproduce la decisión SAA, ese es un escenario
óptimo empírico. El paso siguiente es aprender $q$ (o un sesgo) en función del contexto,
lo que corresponde al pronóstico orientado a la decisión del artículo.

---

## 8. Plan sistemático

| Etapa | Objetivo | Entregable | Criterio de cierre |
|---|---|---|---|
| **F0** Ecuaciones | Fijar el modelo de §3 | Este documento, más la resolución de los puntos de §3.1 | Cada término trazado a un bloque del `.slx` |
| **F1** Réplica algebraica | Escribir el ROM en un lenguaje con diferenciación automática | Código del ROM (16 estados) y su estado estacionario (N0) | La réplica reproduce la PCT de Simulink en el horizonte, con error ≪ 1 K, sobre trayectorias archivadas |
| **F2** Monolítico no lineal (N3) | Resolver P0 como NLP | Solver sobre los estados $t=1161$ y $1180$ s | Acción y costo comparados con Luus–Jaakola, re-simulando en Simulink (V3) |
| **F3** Linealizado (N1, N2) | QP y funciones de transferencia | $A,B,E,C,D$ en cada consulta; $G(s)$ de $\dot m_{\rm ACS}$ a PCT; MPC lineal | Error de N1/N2 frente a N3 ≪ 1 K; brecha de costo medida |
| **F4** Riesgo | Las capas de §7 | Restricción probabilística, SAA y experimento de cuantiles | Una decisión que cambia por riesgo, con intervalos de confianza |
| **F5** Lazo cerrado | Tiempo real y extensiones | 61 consultas con la mejor variante; luego multietapa y dos protecciones | Menos de 120 s por consulta; frontera costo–riesgo fuera de muestra |

**Para la próxima reunión (3–4 semanas)**: F0 cerrado, F1 con la réplica validada en al
menos una trayectoria, y un primer resultado de F2 o F3 en el estado $t=1180$ s.

---

## 9. Decisiones abiertas

1. **ROM o HRM como restricción**: el ROM (16 estados) es suficiente para empezar; el HRM
   de 12 nodos multiplica los estados TH por 12 y elimina la corrección $p_4$.
2. **Herramienta**: JuMP + Ipopt (Julia), CasADi (Python/MATLAB) o fmincon. Ver §6.
3. **PJA y pared cuasiestática**: reducen la rigidez, pero hay que medir su error.
4. **Corrección de la PCT**: mantener $p_4$, reconectar la BiLSTM o ambas.
5. **Norma del seguimiento**: cuadrática (P0) o $L_1$ (compatible con el teorema de §7.1).
