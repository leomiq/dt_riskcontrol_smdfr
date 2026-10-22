# M1: cota de costo y búsqueda por cajas

Fecha: 22 de octubre de 2026. Se mantiene el costo vigente de `fitness_cost.m`; estos resultados no certifican aún el nuevo objetivo estocástico.

## Problema y criterio

Para la consulta fija, sea `U` el costo de un control factible, `L_b` la cota inferior de cada caja de controles y topología, y `L=min_b L_b`. Si todas las cajas cubren el dominio y el oráculo encierra **todas** las trayectorias de cada caja, entonces `L <= J* <= U`. El procedimiento divide la caja con menor `L_b`, evalúa controles factibles y detiene cuando `U-L <= tol` o se agota el presupuesto. La detención no convierte por sí sola un oráculo aproximado en certificado.

`matlab/m1_branch_and_bound.m` implementa esta búsqueda condicional. `bound_fn(box)` debe devolver una cota válida para todos los controles y todas las topologías DET alcanzables dentro de la caja; `cost_fn(u)` debe usar el mismo modelo y objetivo, y declarar factibilidad. El código guarda cajas pendientes, número de divisiones y evaluaciones. Rechaza valores de cota o costo inválidos. La validez física de `bound_fn` y el redondeo dirigido son obligaciones externas pendientes.

Para una topología fija, `certified_cost_lower_bound_box.m` usa intervalos `T_j in [L_j,H_j]` y `m_k in [l_k,h_k]`:

- Seguimiento: distancia al cuadrado de `T_ref` a cada intervalo térmico.
- Caudal: distancia al cuadrado de `[l_(k+1)-h_k, h_(k+1)-l_k]` a `[-5,5]`, multiplicada por `1e-4`.
- Seguridad: bisagra cuadrática evaluada en `L_j`, con límite absoluto o DSM según el modo.

Las cotas de cada topología se combinan tomando el **mínimo**, porque cualquier topología alcanzable puede contener el óptimo. La unión de intervalos no sustituye la cobertura de las secuencias de eventos ni resuelve por sí misma diferencias de longitud de trayectoria.

## Demostración híbrida controlada

El ejemplo de prueba tiene `u in [0,1]` y dos modos: para `u<=0.5`, `T=[u,1-u]`; para `u>0.5`, `T=[u,1.2-u]`. Las imágenes intervalares de ambas expresiones se calculan exactamente en cada intersección de la caja con su región de modo. El óptimo analítico es `J*=0.5` en `u=0.5`: en el primer modo `u^2+(1-u)^2 >= 0.5`; en el segundo `u^2+(1.2-u)^2 >= 0.72`.

`test_m1_branch_and_bound`, con tolerancia absoluta `1e-3`, obtuvo `L=0.499024391`, `U=0.500000000`, brecha `0.000975609`, 47 divisiones y 95 evaluaciones. También rechazó un oráculo que devolvía una cota negativa. Esta es una prueba del procedimiento con encierros analíticos del ejemplo, **no** una certificación del ROM/DET real. Los números de máquina no usan aritmética intervalar con redondeo dirigido.

## Auditoría del piloto real

`m1_pilot_audit` carga sólo los límites y el objetivo del MAT de consulta, más el resultado local R2024b. Reconstruye cada término del costo desde la trayectoria y el control guardados; verifica límites y la factibilidad muestreada.

| Término | Valor |
|---|---:|
| Seguimiento `J_T` | 0.348406945669 |
| Saltos de caudal `J_m` | 0.013619170139 |
| Seguridad `J_s` | 0 |
| Candidato factible `U` | 0.362026115808 |
| Cota global derivada sólo del caudal `L_m` | 0 |
| Brecha disponible `U-L` | 0.362026115808 |

La transformación `CurrentToMassFlow` es monótona no decreciente en corriente dentro del dominio usado: antes del truncamiento, `dQ/dI=C/sqrt(B^2-4A(DeltaP-CI))>0` cuando el discriminante es positivo, y el truncamiento a cero preserva monotonía. Sus valores en los extremos de cada intervalo de corriente encierran el caudal en aritmética exacta. Esos intervalos admiten secuencias de caudal con diferencias menores que 5 kg/s, por lo que `L_m=0`. No hay aún una cota térmica no trivial para **todos** los controles y ramas. El límite global disponible para el piloto sigue siendo `L=0`; la brecha no demuestra cercanía al óptimo. El resultado del piloto sigue discrepando +11.0 % de la referencia R2024a archivada.

## Estado y siguiente trabajo

M1 tiene ahora el algoritmo, un caso híbrido con cobertura demostrada y la auditoría de un candidato real. **No está cerrado para el ROM/DET**. Para cerrarlo hace falta:

1. Congelar la misma versión de modelo, estado inicial, horizonte y criterio de factibilidad en el oráculo de cota y en la evaluación de candidatos.
2. Encerrar el ROM por cajas de control con errores validados, incluyendo sensibilidad a condiciones iniciales y señales exógenas.
3. Enumerar o encerrar cada transición DET alcanzable por caja, con tiempos y longitudes de trayectoria compatibles con el costo vigente.
4. Usar aritmética intervalar o cotas de error de redondeo para que `L_b` sea fiable numéricamente.
5. Reconciliar la discrepancia del piloto R2024b frente a R2024a antes de comparar algoritmos.

## Reproducción

Desde la raíz del repositorio:

```text
matlab -batch "addpath('matlab'); test_certified_cost_lower_bound_box; test_m1_branch_and_bound"
matlab -batch "addpath('matlab'); m1_pilot_audit"
```

El segundo comando escribe `results/m1_pilot_audit.mat`; lee `models/LJ_test_sampled_vel_ACS_2_loop_1180.mat` y `results/pilot_r2024b_20261022.mat`. No ejecuta Simulink.
