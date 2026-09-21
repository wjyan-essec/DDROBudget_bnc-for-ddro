@echo off
setlocal EnableExtensions

rem Usage: test_n10_windows.bat [replicate 1-10 or all] [time limit in seconds]
rem Examples: test_n10_windows.bat 1 60
rem           test_n10_windows.bat all 120
rem Set PYTHON_EXE beforehand to override the Python interpreter.

set "ROOT=%~dp0"
set "INPUT_DIR=%ROOT%instances\originalKP\N10_rc_0.7_rGamma_1_dc_0.035_dGamma_0.0005"
set "STEM=N10_rc_0.7_rGamma_1_dc_0.035_dGamma_0.0005_r_10_100_"
set "REPLICATE=%~1"
if "%REPLICATE%"=="" set "REPLICATE=1"
set "TIME_LIMIT=%~2"
if "%TIME_LIMIT%"=="" set "TIME_LIMIT=60"

if not defined PYTHON_EXE (
    if exist "%ROOT%.venv\Scripts\python.exe" (
        set "PYTHON_EXE=%ROOT%.venv\Scripts\python.exe"
    ) else (
        set "PYTHON_EXE=python"
    )
)

if not exist "%INPUT_DIR%\%STEM%1.txt" (
    echo N10 input files were not found in "%INPUT_DIR%".
    exit /b 1
)

if /I not "%REPLICATE%"=="all" (
    if not exist "%INPUT_DIR%\%STEM%%REPLICATE%.txt" (
        echo Invalid replicate "%REPLICATE%". Use 1-10 or all.
        exit /b 1
    )
)

set "PYTHONPATH=%ROOT%src;%PYTHONPATH%"
"%PYTHON_EXE%" -c "import cplex, gurobipy, numpy; from python_to_cplex_c_api import cplex_c_api_wrapper"
if errorlevel 1 (
    echo Python dependencies are unavailable. Set PYTHON_EXE to the configured Windows Python executable.
    exit /b 1
)

set "RUN_DIR=%ROOT%test_runs\N10_%RANDOM%_%RANDOM%"
if exist "%RUN_DIR%" (
    echo Test directory already exists: "%RUN_DIR%".
    exit /b 1
)
mkdir "%RUN_DIR%"
if errorlevel 1 exit /b 1

rem N10 SOURCE values end in .txt, while the converter expects r_10_100_k.
rem Normalize copies in the test directory; never change the original inputs.
"%PYTHON_EXE%" -c "from pathlib import Path; import re; src=Path(r'%INPUT_DIR%'); dst=Path(r'%RUN_DIR%\normalized'); dst.mkdir(); [(dst / p.name).write_text(re.sub(r'(?m)^(SOURCE[ \t]+\S+?)\.txt[ \t]*$', r'\1', p.read_text(encoding='utf-8')), encoding='utf-8') for p in sorted(src.glob('*.txt'))]"
if errorlevel 1 (
    echo Failed to prepare N10 input copies. Test directory: "%RUN_DIR%".
    exit /b 1
)

"%PYTHON_EXE%" "%ROOT%src\convert_kp_to_bobilib.py" --input-root "%RUN_DIR%\normalized" --output-dir "%RUN_DIR%\converted" --sizes 10 --densities 100 --expected-count 10 --scale 100 >"%RUN_DIR%\conversion.log" 2>&1
if errorlevel 1 (
    echo Conversion failed. See "%RUN_DIR%\conversion.log".
    exit /b 1
)

echo N10 test directory: "%RUN_DIR%"
echo Time limit per solve: %TIME_LIMIT% seconds
pushd "%RUN_DIR%"
if errorlevel 1 exit /b 1

if /I "%REPLICATE%"=="all" (
    for /L %%N in (1,1,10) do (
        call :solve_one %%N
        if errorlevel 1 goto :failed
    )
) else (
    call :solve_one %REPLICATE%
    if errorlevel 1 goto :failed
)

popd
"%PYTHON_EXE%" -c "import csv,re; from pathlib import Path; ref=Path(r'%ROOT%instances\originalKP\N10_solutions.txt').read_text(encoding='utf-8'); expected=dict((name,float(value)) for name,value in re.findall(r'instance=(\S+).*?Objective value = ([0-9.]+)',ref,re.S)); rows=list(csv.DictReader(Path(r'%RUN_DIR%\summary.csv').open(newline='',encoding='utf-8'))); check=lambda row: 'NO_INCUMBENT' if row['has_incumbent']!='1' else 'INVALID_HIGH' if float(row['incumbent'])>expected[row['instance']]+1e-4 else 'MATCH' if abs(float(row['incumbent'])-expected[row['instance']])<=1e-4 else 'BELOW_REFERENCE'; [print(row['instance'],row['method'],'reference='+str(expected[row['instance']]),'incumbent='+row['incumbent'],'bound='+row['bound'],'status='+row['status'],'check='+check(row)) for row in rows]" >"%RUN_DIR%\comparison.txt"
if errorlevel 1 (
    echo Reference comparison failed. Check "%RUN_DIR%\comparison.txt".
    exit /b 1
)
type "%RUN_DIR%\comparison.txt"
echo Finished. Results: "%RUN_DIR%\summary.csv"
echo Comparison: "%RUN_DIR%\comparison.txt"
exit /b 0

:solve_one
set "INSTANCE=%RUN_DIR%\converted\%STEM%%~1.mps"
for %%M in (branchandbound intersection) do (
    echo Running N10 replicate %~1 with %%M...
    "%PYTHON_EXE%" "%ROOT%src\run_bnc.py" --instance_file "%INSTANCE%" --instance_type bobilib --lower_level general --cuts %%M --projected 1 --separation integer --time_lim %TIME_LIMIT% --verbose_level 0 --output_csv "%RUN_DIR%\summary.csv" >"%RUN_DIR%\%STEM%%~1_%%M.log" 2>&1
    if errorlevel 1 (
        echo Solve failed. See "%RUN_DIR%\%STEM%%~1_%%M.log".
        exit /b 1
    )
)
exit /b 0

:failed
popd
echo Test stopped after an error. Completed rows, if any, are in "%RUN_DIR%\summary.csv".
exit /b 1
