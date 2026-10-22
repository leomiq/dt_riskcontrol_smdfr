# Reproducción del piloto y procedencia del ajuste poly4

Fecha: 22 de octubre de 2026. Base original: commit 3266149 de
leomiq/dt_riskcontrol_smdfr. Entorno local: MATLAB R2024b en Windows;
la distribución original documenta MATLAB R2024a.

## Bloqueo encontrado

El primer intento de ejecutar pilot_single_run se detuvo antes de la
optimización. El workspace del modelo ROM cargó models/fit_poly.mat, que
contiene un objeto cfit de Curve Fitting Toolbox. Esa toolbox no está
instalada en esta máquina; MATLAB cargó el objeto como uint32 y falló
coeffvalues(fit_poly), dejando poly_vector sin definir.

## Corrección reproducible

Se mantuvo models/fit_poly.mat intacto. Su SHA-256 es
CDA7E1973A7529AAFBE5B9CF453C4F0B5C258BF3868B3223CA3CE4B34E934C6D.
La propiedad MCOS coeffValues del poly4 contiene, en orden p1 a p5:

    -4.636167450115316e-14
     8.963631660311285e-10
    -6.8727883564717865e-06
     0.025355073868340343
    37.840956461339154

matlab/fit_poly_coefficients.m devuelve esos valores y el workspace ROM
los usa antes de aplicar el desplazamiento original de -2.5 K al
término constante. La verificación independiente
python/verify_fit_poly_coefficients.py lee los cinco doubles del MAT
archivado y los compara exactamente con el código MATLAB. Resultado:
PASS. Esta sustitución elimina la necesidad de deserializar cfit, pero
no modifica el modelo ni reestima el polinomio.

## Piloto de estado fijo

Comando:

    matlab -batch "addpath('matlab'); addpath('models'); pilot_single_run"

Referencia v1.0.0 en README: f_opt = 0.326151 a t_q = 1180 s
(0.326139 en el árbol de referencia). Ejecución local R2024b,
semilla 1, registrada en results/pilot_r2024b_20261022.mat:

| Indicador | Resultado local |
|---|---:|
| f_opt | 0.36202611580770533 |
| Diferencia absoluta frente a 0.326151 | +0.03587511580770531 |
| Diferencia relativa | +10.9995 % |
| Primer control ACS | 56.74999989742171 A |
| Temperatura máxima pronosticada | 1179.5789463721424 K |
| Margen mínimo frente a 1244 K | 64.42105362785765 K |
| Pases LJ | 6, parada temprana |
| Tiempo del optimizador | 828.309508 s |
| Tiempo total de pared reportado | 873.10 s |

La ejecución completa funciona después de la reparación cfit, pero
**el valor numérico v1.0.0 no quedó reproducido en R2024b**. No se
atribuye todavía la diferencia al optimizador: primero hay que comparar
estados iniciales regenerados, salidas ROM/HRM, versiones y MEX con
R2024a, y evaluar el mismo vector de control en ambos entornos.
La búsqueda posterior de cotas debe usar un costo de referencia
congelado y una versión de modelo explícita.

El vector de temperatura guardado tiene rango
[1179.274324533973, 1179.5789463721424] K. Con
$T_{ref}=1179.2712$ K, el mayor valor absoluto de DSM es
0.004754396376 < 1: el candidato satisface la prueba de factibilidad
muestreada de fitness_cost. Para el problema de costo **actual**,
esto da una cota superior U=0.36202611580770533. Sin encierros
certificados sobre todos los controles, la única cota inferior global
disponible es L=0 y el gap U-L=0.36202611580770533; no es un
certificado útil de cercanía al óptimo.

## Verificaciones realizadas

    python python/verify_fit_poly_coefficients.py
    matlab -batch "addpath('matlab'); test_certified_cost_lower_bound_box"

Ambas pruebas pasan. La segunda valida la cota algebraica en cajas
puntuales y envolventes sintéticas; no valida el ROM/DET.
