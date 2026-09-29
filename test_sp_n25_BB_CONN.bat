@echo off
setlocal EnableExtensions EnableDelayedExpansion

rem Usage: test_sp_n25_BB_CONN.bat [total time limit in seconds]
rem Runs all 20 N25 SP instances with connectivity cuts only.
rem The TSP runner subtracts each instance's four-bound preprocessing time.
set "ROOT=%~dp0"
set "INSTANCE_DIR=%ROOT%instances\sp_aggregated_bobilib"
set "STEM=SP_N25_rc_1_rGamma_1_dc_0.035_dGamma_0.0005_SP_N25_ID"
set "METHOD=branchandbound"
set "EXPECTED_METHOD=BB-CONN"
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
    echo Generate it first with:
    echo   "%PYTHON_EXE%" "%ROOT%src\convert_sp_to_bobilib.py"
    exit /b 1
)

for /L %%R in (1,1,20) do (
    set "ID=0%%R"
    set "ID=!ID:~-2!"
    if not exist "%INSTANCE_DIR%\%STEM%!ID!.mps" (
        echo Missing MPS: "%INSTANCE_DIR%\%STEM%!ID!.mps"
        exit /b 1
    )
    if not exist "%INSTANCE_DIR%\%STEM%!ID!.aux" (
        echo Missing AUX: "%INSTANCE_DIR%\%STEM%!ID!.aux"
        exit /b 1
    )
)

set "PYTHONPATH=%ROOT%src;%PYTHONPATH%"
"%PYTHON_EXE%" -c "import cplex,gurobipy,numpy; from python_to_cplex_c_api import cplex_c_api_wrapper as w; assert hasattr(w,'CPXgetcallbacknodex')"
if errorlevel 1 (
    echo Python dependencies or the rebuilt CPLEX callback bridge are unavailable.
    echo Set PYTHON_EXE to the configured Windows Python executable and rebuild in src\python_to_cplex_c_api.
    exit /b 1
)

set "RUN_DIR=%ROOT%test_runs\SP_N25_BB_CONN_%RANDOM%_%RANDOM%"
if exist "%RUN_DIR%" (
    echo Test directory already exists: "%RUN_DIR%"
    exit /b 1
)
mkdir "%RUN_DIR%"
if errorlevel 1 exit /b 1

echo SP N25 aggregated instances: 20
echo Method: %EXPECTED_METHOD%
echo Total time limit per instance: %TIME_LIMIT% seconds
echo Results: "%RUN_DIR%\summary.csv"

set /a RUN_COUNT=0
set "FAILED=0"
for /L %%R in (1,1,20) do (
    call :solve_one %%R
    if errorlevel 1 set "FAILED=1"
)

if not exist "%RUN_DIR%\summary.csv" (
    echo No result CSV was produced. Logs: "%RUN_DIR%"
    exit /b 1
)

"%PYTHON_EXE%" -c "import csv,sys; rows=list(csv.DictReader(open(sys.argv[1],newline='',encoding='utf-8'))); names=[r['instance'] for r in rows]; ok=len(rows)==20 and len(set(names))==20 and all(r['problem']=='TSP' and r['formulation']=='DDRO_BI' and r['linearcons']=='A' and r['method']=='BB-CONN' for r in rows); print('Completed rows:',len(rows),'/20'); sys.exit(0 if ok else 1)" "%RUN_DIR%\summary.csv"
if errorlevel 1 set "FAILED=1"

"%PYTHON_EXE%" -c "import csv,sys; rows=list(csv.DictReader(open(sys.argv[1],newline='',encoding='utf-8'))); refs={'%STEM%01':696.34,'%STEM%02':738.74,'%STEM%03':789.21,'%STEM%04':674.59,'%STEM%05':610.31,'%STEM%06':751.70,'%STEM%07':653.27,'%STEM%08':692.37,'%STEM%09':832.79,'%STEM%10':638.92,'%STEM%11':787.43,'%STEM%12':672.73,'%STEM%13':710.74,'%STEM%14':677.87,'%STEM%15':652.01,'%STEM%16':724.64,'%STEM%17':823.49,'%STEM%18':633.60,'%STEM%19':750.38,'%STEM%20':618.22}; check=lambda r:'NO_INCUMBENT' if r['has_incumbent']=='0' else 'INVALID_TOUR' if not r['tour_valid']=='1' else 'INVALID_LOW' if float(r['incumbent'])<refs[r['instance']]-1e-4 else 'MATCH' if abs(float(r['incumbent'])-refs[r['instance']])<=1e-4 else 'FEASIBLE_ABOVE'; invalid=[r for r in rows if check(r).startswith('INVALID') or (r['status']=='OPTIMAL' and not check(r)=='MATCH')]; [print(r['instance'],'reference='+str(refs[r['instance']]),'incumbent='+r['incumbent'],'bound='+r['bound'],'status='+r['status'],'tour_valid='+r['tour_valid'],'check='+check(r)) for r in rows]; print('Invalid results:',len(invalid)); sys.exit(1 if invalid else 0)" "%RUN_DIR%\summary.csv" >"%RUN_DIR%\comparison.txt"
if errorlevel 1 set "FAILED=1"
type "%RUN_DIR%\comparison.txt"

echo Result CSV: "%RUN_DIR%\summary.csv"
echo Comparison: "%RUN_DIR%\comparison.txt"
echo Logs: "%RUN_DIR%"
if "%FAILED%"=="1" (
    echo At least one run failed, a result row is missing, or a correctness check failed.
    exit /b 1
)
echo All 20 SP N25 BB-CONN runs completed without a correctness violation.
exit /b 0

:solve_one
set /a RUN_COUNT+=1
set "ID=0%~1"
set "ID=!ID:~-2!"
set "INSTANCE=%INSTANCE_DIR%\%STEM%!ID!.mps"
echo.
echo ============================================================
echo Run !RUN_COUNT!/20
echo Instance: %STEM%!ID!
echo Method: %EXPECTED_METHOD%
echo Result CSV: "%RUN_DIR%\summary.csv"
echo ============================================================
"%PYTHON_EXE%" "%ROOT%src\run_bnc_tsp.py" --instance_file "%INSTANCE%" --cuts %METHOD% --projected 1 --separation integer --time_lim %TIME_LIMIT% --cplex_cuts -1 --only_root_node 0 --tolerance 1e-6 --verbose_level 0 --write_lps 0 --output_csv "%RUN_DIR%\summary.csv" >"%RUN_DIR%\%STEM%!ID!_%EXPECTED_METHOD%.log" 2>&1
if errorlevel 1 (
    echo Failed: %STEM%!ID!. See "%RUN_DIR%\%STEM%!ID!_%EXPECTED_METHOD%.log".
    exit /b 1
)
echo Finished: %STEM%!ID!
exit /b 0
