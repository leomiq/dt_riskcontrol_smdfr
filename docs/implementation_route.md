# Ruta de implementación: control estocástico informado por riesgo

Fecha: 22 de octubre de 2026. Repositorio de trabajo: leomiq/dt_riskcontrol_smdfr, rama main, punto de partida 3266149e36787c61915d3b0c735ced43ad6e2595 (v1.0.0). Esta carpeta ya era una copia Git local; se conserva y se desarrolla aquí. No se modificó phd_thesis.

## Estado comprobado

- El piloto de consulta fija usa DET_2, luus_jaakola_det y fitness_cost. GA llama la misma función de costo; el ejemplo principal del scheduler tiene las llamadas LJ/GA comentadas y aplica PID en esa variante. Conviene usar el piloto como línea base de optimización y tratar explícitamente los distintos drivers.
- La función fitness_cost suma error cuadrático de temperatura, penalización por saltos de caudal y penalización por temperatura sobre 1244 K, o por DSM si se suministran dsm_thr y lambda_dsm. La factibilidad se devuelve, pero LJ registra candidatos aun cuando feasible es falso; GA filtra los candidatos no factibles. Hay que unificar la política.
- fitness_cost usa sólo forest_data{1}, concatena muestras de nodos y no pondera rutas por probabilidad. La topología y el número de muestras pueden variar según el control. Por ello el costo actual no es una esperanza estocástica bien definida sobre rutas.
- El PSM por IS-MCDET se calcula en DET_scheduler_postprocessing_perc_online.m a partir de rutas terminales y pesos de importancia, después de la búsqueda. El DSM dentro del costo es determinista y no equivale al PSM probabilístico.
- Los ensayos archivados muestran que la penalización DSM sólo cambió el óptimo en 2 de 49 configuraciones a t=1180 s, con umbral bajo la banda alcanzable; en fase aguda el óptimo ya saturaba el ACS. Véase docs/DSM_PENALTY_RESULTS.md. Esto es evidencia del caso ensayado, no una garantía general.
- El árbol contiene el par [demand_ps, ava_ps] del ACS y ramas de fallo/reparación. Con un solo sistema activo, su disponibilidad no introduce una elección de orden de activación entre sistemas. La integración con protección pasiva requiere nuevos estados, guardas y dependencias causales.

## Formulación objetivo propuesta

En cada tiempo de consulta t, una política causal pi decide u_k y activaciones binarias a_k tras observar la historia disponible. El estado híbrido z_k contiene temperaturas, estado del ACS, estados de los sistemas adicionales y del sistema pasivo. La incertidumbre xi reúne demanda, fallo, reparación, perturbación, error del modelo y condiciones iniciales. El modelo produce T_{k+1}=F(T_k,z_k,u_k,a_k,xi_k), con transiciones de z_k condicionadas por demanda, disponibilidad y tiempos de activación.

Objetivo de investigación:

    min_pi E_p[sum_k (w_T(T_k-T_ref)^2 + w_u c_u(u_k,a_k))]
           + lambda_R R_pi
    sujeto a P_p(max_k T_k > U_y | informacion_t) <= epsilon
              y restricciones físicas, lógicas y de no anticipación.

Comparar R_pi = P_p(max T > U_y), CVaR_alpha(max T-U_y), un PSM de cuantil superior y el mínimo DSM por ruta. El evento y la medida de probabilidad deben fijarse antes de comparar métodos. El PSM no debe confundirse con la probabilidad de excedencia: son cantidades relacionadas, pero con unidades e interpretación distintas. En IS, usar pesos p/q de trayectorias completas y registrar tamaño efectivo de muestra, varianza e intervalos de confianza. La distribución de fallo debe ser condicional al estado y a la acción si la acción cambia demanda o carga del componente. Mantener no anticipación: las acciones previas a una bifurcación son idénticas en ramas indistinguibles.

## Hitos y criterios de cierre

### M0. Congelar referencia y reproducibilidad

Registrar SHA, dependencias, semillas, entrada t_q=1180 s, valor de referencia f=0.326151 y tiempos de ejecución. Ejecutar piloto y capturar descomposición de J y secuencia de estados. Separar resultados nuevos de resultados publicados; registrar si no hay licencia MATLAB/Simulink operativa.

Cierre: un comando reproducible, archivo de parámetros y reporte de diferencia respecto al piloto publicado.

### M1. Cota inferior para el costo vigente — iniciado en este commit

Derivar la descomposición J=J_T+J_m+J_s y una cota por cajas. Para cada muestra futura, encerrar T_j en [L_j,U_j]; para cada control, encerrar caudal m_k en [l_k,h_k]. El mínimo de cada término no negativo sobre esas cajas da una cota inferior. El módulo matlab/certified_cost_lower_bound_box.m implementa la fórmula y su prueba sintética verifica cajas puntuales, inclusión y errores de entrada. La cota es **condicional**: sólo certifica el problema si los intervalos cubren todas las trayectorias de la celda de control y la topología DET correspondiente está fija o sus alternativas están todas enumeradas.

Siguiente paso de M1: construir enclosures verificables del ROM/DET (aritmética intervalar, relajación validada o partición con cotas de error certificadas), particionar el dominio de control por topología, calcular LB por celda y UB con un candidato factible evaluado con el mismo costo. Detenerse con (UB-LB)/max(1,abs(UB)) <= tolerancia, conservando informe de celdas pendientes. Un barrido denso o muchas semillas LJ/GA son sólo evidencia empírica, nunca prueba de optimalidad global. Si el simulador no permite enclosures certificables, reportar la cota universal J>=0 y un gap empírico claramente rotulado.

Cierre: LB demostrada, UB factible, gap certificado o constancia explícita de por qué no se obtuvo, más sensibilidad a semillas y presupuesto computacional. No comparar costos de topologías con distinta cantidad de muestras sin una convención de horizonte común.

### M2. Evaluación de riesgo dentro de cada candidato

Extraer el postprocesamiento PSM a una función pura que recibe rutas, probabilidades objetivo/propuesta y temperatura futura. Verificar normalización p/q, rutas sin fallos, fallos raros, pesos extremos y horizonte futuro. Envolver DET_2 para devolver una tabla de escenarios con ID, secuencia de eventos, probabilidad, T_k y estado del sistema. Hacer que cada evaluación de pi calcule J y R_pi usando las **mismas** realizaciones de incertidumbre (common random numbers); usar muestreo adaptativo y confianza simultánea si se comparan muchos candidatos. Primero comparar penalización de riesgo, restricción de probabilidad y CVaR; después escoger por frontera costo–riesgo y coste de cómputo.

Cierre: cambio de decisión atribuible a riesgo calculado para cada candidato, resultados con intervalos de confianza y reproducibilidad. Una penalización DSM que no usa probabilidades no cumple este hito.

### M3. Optimización estocástica y validación

Probar SAA con escenarios fijos, búsqueda local/multistart, optimización por escenarios y, si el tamaño lo permite, programación entera o descomposición para decisiones binarias. Evaluar costo fuera de muestra y riesgo con semillas independientes; calibrar surrogate antes de usarlo para búsqueda; medir latencia frente al periodo de control. Comparar contra LJ, GA, PID y la política DSM bajo igual presupuesto.

Cierre: frontera de Pareto fuera de muestra, incertidumbre de estimación y caso de falla documentado.

### M4. Varios sistemas de protección y seguridad pasiva

Añadir al menos dos sistemas activos con estados disponible, demandado, activado, fallado y reparado; representar retardos, incompatibilidades, recursos compartidos y prioridad. El sistema pasivo tiene demanda y disponibilidad propias, con dependencia condicional de la secuencia activa. Generar árbol de eventos que preserve orden temporal y estados absorbentes; podar sólo con regla auditable y cota del riesgo omitido. La política debe elegir activación y control continuo condicionados a observaciones, sin conocer fallos futuros.

Cierre: casos sintéticos que distingan A->B, B->A, fallos comunes, ACS fallido y activación pasiva; comparación de políticas y trazas causales.

## Decisiones abiertas

1. Definir la métrica principal: probabilidad de excedencia, PSM de cuantil o CVaR, nivel alfa/epsilon y tiempo de misión.
2. Elegir sistemas adicionales, sus tasas/leyes de fallo, retardos y dependencias con el sistema pasivo.
3. Fijar tolerancia del gap, presupuesto por consulta y criterio de incertidumbre estadística.
4. Acordar si el costo debe ser esperanza por ruta, peor caso, o riesgo más costo esperado. La suma actual de nodos sin ponderar es una línea base de software, no el objetivo final.

## Comandos de verificación

    matlab -batch "addpath('matlab'); test_certified_cost_lower_bound_box"
    matlab -batch "addpath('matlab'); addpath('models'); pilot_single_run"

El segundo comando requiere MATLAB R2024a, Simulink, modelos y MEX; escribe pilot_result.mat en el directorio de ejecución. Ejecutarlo desde un directorio de resultados dedicado.
