# Checks that the tools this project needs are installed and reachable.
# Run from the repository root:  powershell -ExecutionPolicy Bypass -File tools\check_env.ps1
$ok = $true
function Report($name, $found, $detail) {
    $mark = if ($found) { 'OK  ' } else { 'FALTA' }
    '{0}  {1,-28} {2}' -f $mark, $name, $detail
    if (-not $found) { $script:ok = $false }
}

$git = Get-Command git -ErrorAction SilentlyContinue
Report 'Git' ($null -ne $git) $(if ($git) { (git --version) } else { 'winget install Git.Git' })
if ($git) {
    $u = git config user.name
    Report 'git user.name / email' ([bool]$u) $(if ($u) { "$u <$(git config user.email)>" } else { 'git config --global user.name ...' })
}

$ml = Get-Command matlab -ErrorAction SilentlyContinue
Report 'MATLAB en PATH' ($null -ne $ml) $(if ($ml) { $ml.Source } else { 'instalar R2024b y agregar ...\MATLAB\R2024b\bin al PATH' })

$jl = Get-Command julia -ErrorAction SilentlyContinue
$jv = if ($jl) { (julia --version) } else { 'winget install Julialang.Juliaup; juliaup add 1.10.4' }
Report 'Julia' ($null -ne $jl) $jv
if ($jl) {
    $inst = Test-Path (Join-Path $env:USERPROFILE '.julia\packages\JuMP')
    Report 'Paquetes Julia (JuMP...)' $inst $(if ($inst) { 'instalados' } else { "julia --project=julia -e 'using Pkg; Pkg.instantiate()'" })
}

$venv = Test-Path '.venv\Scripts\python.exe'
$pv = if ($venv) { (& .venv\Scripts\python.exe --version) } else { 'python -m venv .venv (Python 3.11)' }
Report 'Python venv (.venv)' $venv $pv
if ($venv) {
    foreach ($m in 'numpy', 'scipy', 'matplotlib', 'h5py', 'torch') {
        & .venv\Scripts\python.exe -c "import $m" 2>$null
        $has = ($LASTEXITCODE -eq 0)
        $tip = if ($m -eq 'torch') { 'opcional: requirements-bilstm.txt' } else { 'pip install -r python/requirements.txt' }
        Report "  python: $m" $has $(if ($has) { '' } else { $tip })
    }
}

$gpu = Get-Command nvidia-smi -ErrorAction SilentlyContinue
$gd = if ($gpu) { (nvidia-smi --query-gpu=name,driver_version --format=csv,noheader) } else { 'sin GPU NVIDIA: torch en CPU' }
Report 'GPU NVIDIA (opcional)' $true $gd

foreach ($f in 'models\LJ_test_sampled_vel_ACS_2_loop_1180.mat', 'models\rom_MCDET_v2_2022b_control_CMEX.slx',
               'models\SMDFR_hrm_MCDET_2022b_control_CMEX.slx', 'matlab\msfcn_schedule_rom_v2.mexw64') {
    Report "archivo $f" (Test-Path $f) ''
}
$pth = Test-Path 'models\model_100x100_drop2_BILSTM.pth'
'{0}  {1,-28} {2}' -f $(if ($pth) { 'OK  ' } else { 'AVISO' }), 'pesos BiLSTM (mc_inference)', $(if ($pth) { '' } else { 'model_100x100_drop2_BILSTM.pth no esta en el repo (ver docs/setup_entorno.md)' })

$ram = [math]::Round((Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory / 1GB, 1)
'INFO  RAM                          {0} GB' -f $ram
if ($ok) { "`nEntorno completo." } else { "`nFaltan componentes: ver docs/setup_entorno.md" }
