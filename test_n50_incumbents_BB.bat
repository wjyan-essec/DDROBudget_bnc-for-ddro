@echo off
setlocal EnableExtensions EnableDelayedExpansion

rem Usage: test_n50_incumbents_BB.bat [total time limit in seconds]
rem Re-tests with BB the five N50 instances that produced incumbents previously.
set "ROOT=%~dp0"
set "INSTANCE_DIR=%ROOT%instances\kp_disaggregated_bobilib"
set "STEM=N50_rc_0.7_rGamma_1_dc_0.035_dGamma_0.0005_jeu_50_"
set "METHOD=branchandbound"
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

call :check_instance 50 3
if errorlevel 1 exit /b 1
call :check_instance 50 4
if errorlevel 1 exit /b 1
call :check_instance 75 1
if errorlevel 1 exit /b 1
call :check_instance 75 2
if errorlevel 1 exit /b 1
call :check_instance 75 5
if errorlevel 1 exit /b 1

set "PYTHONPATH=%ROOT%src;%PYTHONPATH%"
"%PYTHON_EXE%" -c "import cplex, gurobipy, numpy; from python_to_cplex_c_api import cplex_c_api_wrapper"
if errorlevel 1 (
    echo Python dependencies are unavailable. Set PYTHON_EXE to the configured Windows Python executable.
    exit /b 1
)

set "RUN_DIR=%ROOT%test_runs\N50_INCUMBENTS_BB_%RANDOM%_%RANDOM%"
if exist "%RUN_DIR%" (
    echo Test directory already exists: "%RUN_DIR%"
    exit /b 1
)
mkdir "%RUN_DIR%"
if errorlevel 1 exit /b 1

echo N50 selected incumbent instances: 5
echo Method: %METHOD%
echo Total time limit per instance: %TIME_LIMIT% seconds
echo Results: "%RUN_DIR%\summary.csv"

set /a RUN_COUNT=0
set "FAILED=0"
call :solve_one 50 3
if errorlevel 1 set "FAILED=1"
call :solve_one 50 4
if errorlevel 1 set "FAILED=1"
call :solve_one 75 1
if errorlevel 1 set "FAILED=1"
call :solve_one 75 2
if errorlevel 1 set "FAILED=1"
call :solve_one 75 5
if errorlevel 1 set "FAILED=1"

if not exist "%RUN_DIR%\summary.csv" (
    echo No result CSV was produced. Logs: "%RUN_DIR%"
    exit /b 1
)

"%PYTHON_EXE%" -c "import csv,sys; rows=list(csv.DictReader(open(sys.argv[1],newline='',encoding='utf-8'))); names=[r['instance'] for r in rows]; print('Completed rows:',len(rows),'/5'); sys.exit(0 if len(rows)==5 and len(set(names))==5 and all(r['method']=='branchandbound' for r in rows) else 1)" "%RUN_DIR%\summary.csv"
if errorlevel 1 set "FAILED=1"

"%PYTHON_EXE%" -c "import csv,sys; rows=list(csv.DictReader(open(sys.argv[1],newline='',encoding='utf-8'))); refs={'%STEM%50_3':197.28,'%STEM%50_4':353.41,'%STEM%75_1':126.78,'%STEM%75_2':150.94,'%STEM%75_5':140.22}; bad=[r for r in rows if r['has_incumbent']=='1' and float(r['incumbent'])>refs[r['instance']]+1e-6]; [print(r['instance'],'incumbent='+r['incumbent'],'reference='+str(refs[r['instance']]),'INVALID_HIGH' if r in bad else 'OK') for r in rows if r['has_incumbent']=='1']; print('Invalid high incumbents:',len(bad)); sys.exit(1 if bad else 0)" "%RUN_DIR%\summary.csv" >"%RUN_DIR%\comparison.txt"
if errorlevel 1 set "FAILED=1"
type "%RUN_DIR%\comparison.txt"

echo Result CSV: "%RUN_DIR%\summary.csv"
echo Comparison: "%RUN_DIR%\comparison.txt"
echo Logs: "%RUN_DIR%"
if "%FAILED%"=="1" (
    echo At least one run failed, a result row is missing, or an invalid high incumbent was found.
    exit /b 1
)
echo All 5 selected N50 BB runs completed without an invalid high incumbent.
exit /b 0

:check_instance
if not exist "%INSTANCE_DIR%\%STEM%%~1_%~2.mps" (
    echo Missing MPS: "%INSTANCE_DIR%\%STEM%%~1_%~2.mps"
    exit /b 1
)
if not exist "%INSTANCE_DIR%\%STEM%%~1_%~2.aux" (
    echo Missing AUX: "%INSTANCE_DIR%\%STEM%%~1_%~2.aux"
    exit /b 1
)
exit /b 0

:solve_one
set /a RUN_COUNT+=1
set "INSTANCE=%INSTANCE_DIR%\%STEM%%~1_%~2.mps"
echo.
echo ============================================================
echo Run !RUN_COUNT!/5
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
