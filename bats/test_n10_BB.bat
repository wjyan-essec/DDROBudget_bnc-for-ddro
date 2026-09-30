@echo off
setlocal EnableExtensions

rem Usage: test_n10_bb_all.bat [time limit in seconds]
rem Example: test_n10_bb_all.bat 60
rem Set PYTHON_EXE beforehand to override the Python interpreter.

set "ROOT=%~dp0"
set "INSTANCE_DIR=%ROOT%instances\kp_disaggregated_bobilib_n10"
set "STEM=N10_rc_0.7_rGamma_1_dc_0.035_dGamma_0.0005_r_10_100_"
set "TIME_LIMIT=%~1"
if "%TIME_LIMIT%"=="" set "TIME_LIMIT=300"

if not defined PYTHON_EXE (
    if exist "%ROOT%.venv\Scripts\python.exe" (
        set "PYTHON_EXE=%ROOT%.venv\Scripts\python.exe"
    ) else (
        set "PYTHON_EXE=python"
    )
)

if not exist "%INSTANCE_DIR%\conversion_manifest.csv" (
    echo N10 conversion manifest was not found in "%INSTANCE_DIR%".
    exit /b 1
)

set "PYTHONPATH=%ROOT%src;%PYTHONPATH%"
"%PYTHON_EXE%" -c "import cplex, gurobipy, numpy; from python_to_cplex_c_api import cplex_c_api_wrapper"
if errorlevel 1 (
    echo Python dependencies are unavailable. Set PYTHON_EXE to the configured Windows Python executable.
    exit /b 1
)

set "RUN_DIR=%ROOT%test_runs\N10_BB_%RANDOM%_%RANDOM%"
if exist "%RUN_DIR%" (
    echo Test directory already exists: "%RUN_DIR%".
    exit /b 1
)
mkdir "%RUN_DIR%"
if errorlevel 1 exit /b 1

echo N10 instances: "%INSTANCE_DIR%"
echo Time limit per solve: %TIME_LIMIT% seconds
echo Results: "%RUN_DIR%\summary.csv"

set "FAILED=0"
for /L %%N in (1,1,10) do (
    call :solve_one %%N
    if errorlevel 1 set "FAILED=1"
)

if not exist "%RUN_DIR%\summary.csv" (
    echo No result CSV was produced. Logs: "%RUN_DIR%"
    exit /b 1
)

"%PYTHON_EXE%" -c "import csv,re,sys; from pathlib import Path; ref=Path(r'%ROOT%instances\originalKP\N10_solutions.txt').read_text(encoding='utf-8'); expected=dict((name,float(value)) for name,value in re.findall(r'instance=(\S+).*?Objective value = ([0-9.]+)',ref,re.S)); rows=list(csv.DictReader(Path(r'%RUN_DIR%\summary.csv').open(newline='',encoding='utf-8'))); check=lambda row: 'NO_INCUMBENT' if row['has_incumbent']!='1' else 'INVALID_HIGH' if float(row['incumbent'])>expected[row['instance']]+1e-4 else 'MATCH' if abs(float(row['incumbent'])-expected[row['instance']])<=1e-4 else 'BELOW_REFERENCE'; [print(row['instance'],'reference='+str(expected[row['instance']]),'incumbent='+row['incumbent'],'bound='+row['bound'],'status='+row['status'],'check='+check(row)) for row in rows]; print('completed_rows='+str(len(rows))+'/10'); sys.exit(0 if len(rows)==10 and all(check(row)=='MATCH' for row in rows) else 1)" >"%RUN_DIR%\comparison.txt"
if errorlevel 1 set "FAILED=1"

type "%RUN_DIR%\comparison.txt"
echo Result CSV: "%RUN_DIR%\summary.csv"
echo Comparison: "%RUN_DIR%\comparison.txt"
echo Logs: "%RUN_DIR%"

if "%FAILED%"=="1" (
    echo At least one run failed or did not match the reference.
    exit /b 1
)
echo All 10 BB results match the reference.
exit /b 0

:solve_one
set "INSTANCE=%INSTANCE_DIR%\%STEM%%~1.mps"
if not exist "%INSTANCE%" (
    echo Missing MPS for replicate %~1: "%INSTANCE%"
    exit /b 1
)
if not exist "%INSTANCE_DIR%\%STEM%%~1.aux" (
    echo Missing AUX for replicate %~1.
    exit /b 1
)

echo Running N10 replicate %~1 with branchandbound...
"%PYTHON_EXE%" "%ROOT%src\run_bnc.py" --instance_file "%INSTANCE%" --instance_type bobilib --lower_level general --cuts branchandbound --projected 1 --separation integer --time_lim %TIME_LIMIT% --verbose_level 0 --output_csv "%RUN_DIR%\summary.csv" >"%RUN_DIR%\%STEM%%~1_branchandbound.log" 2>&1
if errorlevel 1 (
    echo Solve failed for replicate %~1. See "%RUN_DIR%\%STEM%%~1_branchandbound.log".
    exit /b 1
)
echo Finished replicate %~1. CSV: "%RUN_DIR%\summary.csv"
exit /b 0