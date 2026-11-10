# Entorno de trabajo: instalación en un PC nuevo

Probado en Windows 11 x64. Las versiones indicadas son las que produjeron los resultados de
noviembre de 2026 (F1–F3, E2).

## Resumen

| Componente | Versión | ¿Para qué? | ¿Obligatorio? |
|---|---|---|---|
| Git for Windows | reciente | control de versiones | sí |
| MATLAB + Simulink | **R2024b** (el archivo de referencia es de R2024a) | HRM, ROM, DET, LJ/GA, lazo cerrado | sí |
| Parallel Computing Toolbox | misma versión | `parfor` de LJ y GA (sin ella corren en serie, ~8× más lento) | recomendado |
| Statistics and Machine Learning Toolbox | misma versión | `prctile` en el lazo completo | sólo `run_phase2c_full_loop` |
| MinGW-w64 (add-on de MATLAB) | — | recompilar las S-functions | no: los `.mexw64` están en el repo |
| Julia | **1.10.4** (fijado en `julia/Manifest.toml`) | réplica del ROM, NLP (JuMP + Ipopt), QP (Clarabel), cotas por intervalos (IntervalArithmetic) | sí, para F1–F4 |
| Python | **3.11** (MATLAB R2024b admite 3.9–3.12; **no 3.13**) | figuras, KPI, inferencia BiLSTM | sí |
| PyTorch con CUDA 12.1 | ≥ 2.2 | MC-dropout de la BiLSTM (en CPU basta para F4) | sólo lazo completo y F4 |
| Tectonic (o XeLaTeX) | 0.15.0 | compilar `docs/estado_mpc_monolitico.tex` con la fuente Inter (`fontspec`) | sólo para regenerar el PDF |

Ya **no** hace falta la Curve Fitting Toolbox (ver `docs/pilot_reproducibility.md`).

## 1. Git

```powershell
winget install --id Git.Git -e
# abrir una terminal nueva
git config --global user.name  "Leonardo Andrés Miqueles Leiva"
git config --global user.email "leonardomiqueles@gmail.com"
git config --global core.autocrlf true
```

**Sobre OneDrive.** El repositorio vive en una carpeta sincronizada por OneDrive, con `.git`
incluido. Funciona, pero OneDrive puede corromper `.git` si dos PC escriben a la vez, y
sincroniza también los cachés de Simulink. Recomendación: clonar en una carpeta local fuera
de OneDrive (p. ej. `C:\dev\dt_riskcontrol_smdfr`) desde el remoto de GitHub, y copiar a
mano lo que no está versionado (`20261011_docs/`, archivos `.mat` sueltos que necesites).

**Ya ocurrió (6-11-2026)**: con las dos PC abiertas sobre la carpeta, OneDrive reemplazó
`README.md`, `docs/monolithic_formulation_plan.md` y `julia/Project.toml`/`Manifest.toml`
por versiones más antiguas de la otra máquina. Guardó las de leo_acer como copias
`*-leo_acer.*`, que hubo que restaurar a mano. Antes de trabajar, revisar `git status` en
busca de archivos `*-<nombre-del-PC>.*` o de cambios que no son propios.

Si se sigue trabajando en la copia de OneDrive, marcar la carpeta del proyecto con
**"Mantener siempre en este dispositivo"** (clic derecho en el Explorador). Con *Archivos a
petición*, los `.mat`, `.slx` y `.git` pueden quedar sólo en la nube y fallar al abrirse
desde MATLAB, Julia o Git ("La operación de nube no es válida"). Si
Git avisa *dubious ownership* al usar la copia de OneDrive:

```powershell
git config --global --add safe.directory "C:/Users/leona/OneDrive/Escritorio/Current Projects/dt_riskcontrol_smdfr"
```

## 2. MATLAB R2024b

1. Instalar MATLAB R2024b con **Simulink**, **Parallel Computing Toolbox** y
   **Statistics and Machine Learning Toolbox**.
2. Agregar `C:\Program Files\MATLAB\R2024b\bin` al `PATH` del usuario, para que funcione
   `matlab -batch` desde la terminal (los scripts y Claude Code lo usan así).
3. Borrar los cachés de compilación que vienen de otro PC o versión (se regeneran solos y
   están en `.gitignore`):

```powershell
Remove-Item -Recurse -Force slprj, *.slxc -ErrorAction SilentlyContinue
```

4. Prueba rápida (≈ 15 min con `parfor` en serie):

```powershell
matlab -batch "addpath('matlab'); addpath('models'); pilot_single_run"
```

Valor esperado en R2024b: `f = 0.362026` (ver `docs/pilot_reproducibility.md`).

## 3. Julia 1.10.4

```powershell
winget install --id Julialang.Juliaup -e
# terminal nueva
juliaup add 1.10.4
juliaup default 1.10.4
julia --project=julia -e 'using Pkg; Pkg.instantiate()'
```

Todos los comandos se ejecutan **desde la carpeta del repositorio** (`cd` a ella primero):
`--project=julia` es una ruta relativa. Sin `juliaup add 1.10.4`, juliaup instala la última
versión estable (1.13 en nov. 2026), que no corresponde al `Manifest.toml`.

`Pkg.instantiate()` instala exactamente las versiones de `julia/Manifest.toml` (JuMP 1.32.1,
Ipopt 1.16, ForwardDiff, MAT). En PowerShell 5.1 no se deben poner comillas dobles dentro
del `-e`: se pierden al pasar a `julia.exe`. Para código más largo, usar un archivo `.jl`.

Prueba (segundos):

```powershell
julia --project=julia julia/validate_replica.jl
julia --project=julia julia/f2_monolithic_nlp.jl 1180
```

Esperado: error de la réplica ≈ 2.2·10⁻⁶ K en el horizonte; NLP con J = 0.244.

## 4. Python 3.11

```powershell
winget install --id Python.Python.3.11 -e
```

Desactivar el alias de la Microsoft Store (Configuración → Aplicaciones → Configuración
avanzada de aplicaciones → Alias de ejecución de aplicaciones → apagar `python.exe` y
`python3.exe`), o llamar siempre a `py -3.11`.

```powershell
py -3.11 -m venv .venv
.venv\Scripts\python -m pip install --upgrade pip
.venv\Scripts\python -m pip install -r python/requirements.txt
```

### 4.1 BiLSTM (sólo para el lazo completo)

```powershell
nvidia-smi    # driver >= 530 para CUDA 12.1
.venv\Scripts\python -m pip install -r python/requirements-bilstm.txt --index-url https://download.pytorch.org/whl/cu121 --extra-index-url https://pypi.org/simple
```

Conectar MATLAB a ese Python (una vez; MATLAB lo recuerda). `OutOfProcess` evita choques de
DLL entre torch y MATLAB:

```matlab
pyenv('Version', fullfile(pwd, '.venv', 'Scripts', 'python.exe'), 'ExecutionMode', 'OutOfProcess');
insert(py.sys.path, int32(0), fullfile(pwd, 'python'));
addpath('matlab'); addpath('models');
check_torch_env
```

**Pendiente conocido**: `python/mc_inference.py` carga `model_100x100_drop2_BILSTM.pth`, que
**no está en el repositorio** (sólo está `models/model_200x200_drop1_1.pth`). Hay que
recuperarlo del PC anterior o de Cineca, o comprobar que la arquitectura de
`model_200x200_drop1_1.pth` coincide y cambiar la ruta. Además, el archivo se busca en el
directorio de trabajo de Python: ejecutar MATLAB desde la carpeta que lo contiene o usar
una ruta absoluta.
Avance (6-11-2026): la arquitectura de `model_200x200_drop1_1.pth` **sí coincide**, y su
ensamble contiene la muestra archivada en 1180 s (ver §8.2 de
`monolithic_formulation_plan.md`). Los scripts F4 (`python/f4_*.py`) lo cargan con ruta
absoluta. En leo_acer (sin GPU) basta torch en CPU:
`.venv\Scripts\python -m pip install "torch>=2.2" --index-url https://download.pytorch.org/whl/cpu`.

### 4.2 Documento de estado (LaTeX)

`docs/estado_mpc_monolitico.tex` usa `fontspec` con la fuente Inter
(`assets/fonts/InterVariable.ttf`), así que necesita XeLaTeX o LuaLaTeX; `pdflatex` no
sirve. En leo_acer se usa **Tectonic 0.15.0**: un ejecutable único basado en XeTeX, en
`C:\Users\leona\bin\tectonic.exe`, que descarga los paquetes de LaTeX la primera vez. Las
figuras se regeneran desde los resultados guardados, sin MATLAB ni Julia:

```powershell
.venv\Scripts\python python/make_status_figures.py     # -> docs/figures/*.pdf
cd docs; tectonic -X compile estado_mpc_monolitico.tex
```

Detalle conocido: con Tectonic, `\setmainfont` no admite `UprightFont=` junto con `Path=`
(no encuentra la fuente), y babel en español choca con `\max` y `\min` dentro de los pies
de figura: hay que usar `\rm max` y `\rm min`.

## 5. Claude Code

- La memoria del proyecto se guarda por usuario y ruta en
  `%USERPROFILE%\.claude\projects\<ruta-del-proyecto>\memory\`. Si cambia la ruta del
  repositorio (p. ej. al clonar fuera de OneDrive), copiar esa carpeta desde el PC anterior.
- Con 16 GB de RAM, Claude Code detuvo varias corridas largas de MATLAB en segundo plano por
  falta de memoria. Para evitarlo: correr un solo MATLAB a la vez, cerrar aplicaciones
  pesadas y, si hace falta, iniciar Claude Code con la variable de entorno
  `CLAUDE_CODE_DISABLE_BG_SHELL_PRESSURE_REAP=1`.

## 6. PowerShell

Si al abrir PowerShell aparece *"no se puede cargar ... profile.ps1 porque la ejecución de
scripts está deshabilitada"*, es sólo el perfil personal (en `OneDrive\Documentos`) que no
se carga; no afecta al proyecto. Para permitir scripts locales (también `tools\check_env.ps1`):

```powershell
Set-ExecutionPolicy -Scope CurrentUser -ExecutionPolicy RemoteSigned
```

## 7. Verificación

```powershell
powershell -ExecutionPolicy Bypass -File tools\check_env.ps1
```

Lista cada componente como `OK` o `FALTA`, con el comando para instalarlo.
