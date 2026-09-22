@echo off
setlocal EnableExtensions EnableDelayedExpansion

rem Usage: test_n50_BC.bat [total time limit in seconds]
rem The runner subtracts each instance's bound time from this total limit.
set "ROOT=%~dp0"
set "INSTANCE_DIR=%ROOT%instances\kp_disaggregated_bobilib"
set "STEM=N50_rc_0.7_rGamma_1_dc_0.035_dGamma_0.0005_jeu_50_"
set "METHOD=intersection"
set "TIME_LIMIT=%~1"
if "%TIME_LIMIT%"=="" set "TIME_LIMIT=3600"

if not defined PYTHON_EXE (
    if exist "%ROOT%.venv\Scripts\python.exe" (
        set "PYTHON_EXE=%ROOT%.venv\Scripts\python.exe"
    ) else (
        set "PYTHON_EXE=python"
    )
)

if not exist "%INSTANCE_DIR%\conversion_manifest.csv" (
    echo Missing conversion manifest: "%INSTANCE_DIR%\conversion_manifest.csv"
    exit /b 1
)

for %%D in (50 75) do (
    for /L %%R in (1,1,10) do (
        if not exist "%INSTANCE_DIR%\%STEM%%%D_%%R.mps" (
            echo Missing MPS: "%INSTANCE_DIR%\%STEM%%%D_%%R.mps"
            exit /b 1
        )
        if not exist "%INSTANCE_DIR%\%STEM%%%D_%%R.aux" (
            echo Missing AUX: "%INSTANCE_DIR%\%STEM%%%D_%%R.aux"
            exit /b 1
        )
    )
)

set "PYTHONPATH=%ROOT%src;%PYTHONPATH%"
"%PYTHON_EXE%" -c "import cplex, gurobipy, numpy; from python_to_cplex_c_api import cplex_c_api_wrapper"
if errorlevel 1 (
    echo Python dependencies are unavailable. Set PYTHON_EXE to the configured Windows Python executable.
    exit /b 1
)

set "RUN_DIR=%ROOT%test_runs\N50_BC_%RANDOM%_%RANDOM%"
if exist "%RUN_DIR%" (
    echo Test directory already exists: "%RUN_DIR%"
    exit /b 1
)
mkdir "%RUN_DIR%"
if errorlevel 1 exit /b 1

echo N50 density 50 and 75: 20 instances
echo Method: %METHOD%
echo Total time limit per instance: %TIME_LIMIT% seconds
echo Results: "%RUN_DIR%\summary.csv"

set /a RUN_COUNT=0
set "FAILED=0"
for %%D in (50 75) do (
    for /L %%R in (1,1,10) do (
        call :solve_one %%D %%R
        if errorlevel 1 set "FAILED=1"
    )
)

if not exist "%RUN_DIR%\summary.csv" (
    echo No result CSV was produced. Logs: "%RUN_DIR%"
    exit /b 1
)

"%PYTHON_EXE%" -c "import csv,sys; rows=list(csv.DictReader(open(sys.argv[1],newline='',encoding='utf-8'))); names=[r['instance'] for r in rows]; print('Completed rows:',len(rows),'/20'); sys.exit(0 if len(rows)==20 and len(set(names))==20 else 1)" "%RUN_DIR%\summary.csv"
if errorlevel 1 set "FAILED=1"

echo Result CSV: "%RUN_DIR%\summary.csv"
echo Logs: "%RUN_DIR%"
if "%FAILED%"=="1" (
    echo At least one run failed or a result row is missing.
    exit /b 1
)
echo All 20 N50 intersection runs completed.
exit /b 0

:solve_one
set /a RUN_COUNT+=1
set "INSTANCE=%INSTANCE_DIR%\%STEM%%~1_%~2.mps"
echo.
echo ============================================================
echo Run !RUN_COUNT!/20
echo Instance: %STEM%%~1_%~2
echo Method: %METHOD%
echo Result CSV: "%RUN_DIR%\summary.csv"
echo ============================================================
"%PYTHON_EXE%" "%ROOT%src\run_bnc.py" --instance_file "%INSTANCE%" --instance_type bobilib --lower_level general --cuts %METHOD% --projected 1 --separation integer --time_lim %TIME_LIMIT% --max_cuts 20 --cplex_cuts -1 --only_root_node 0 --tolerance 1e-6 --verbose_level 0 --write_lps 0 --output_csv "%RUN_DIR%\summary.csv" >"%RUN_DIR%\%STEM%%~1_%~2_%METHOD%.log" 2>&1
if errorlevel 1 (
    echo Failed: %STEM%%~1_%~2. See "%RUN_DIR%\%STEM%%~1_%~2_%METHOD%.log".
    exit /b 1
)
echo Finished: %STEM%%~1_%~2
exit /b 0
