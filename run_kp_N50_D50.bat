@echo off
setlocal EnableExtensions EnableDelayedExpansion

REM Always run from the repository root.
pushd "%~dp0"

set "INSTANCE_DIR=instances\kp_disaggregated_bobilib"
set "TIME_LIMIT=3600"

REM Create a unique CSV for each batch run.
for /f %%I in ('powershell -NoProfile -Command "Get-Date -Format yyyyMMdd_HHmmss"') do set "TIMESTAMP=%%I"
set "OUTPUT_CSV=results\kp_N50_D50_!TIMESTAMP!.csv"

if not exist "results" mkdir "results"

REM Verify that matching instances exist.
dir /b "%INSTANCE_DIR%\N50_*_jeu_50_50_*.mps" >nul 2>&1
if errorlevel 1 (
    echo ERROR: No N50 D50 instances were found in:
    echo %INSTANCE_DIR%
    popd
    exit /b 1
)

set /a RUN_COUNT=0

for %%F in ("%INSTANCE_DIR%\N50_*_jeu_50_50_*.mps") do (
    for %%M in (branchandbound intersection) do (
        set /a RUN_COUNT+=1

        echo.
        echo ============================================================
        echo Run !RUN_COUNT!
        echo Instance: %%~nF
        echo Method:   %%M
        echo ============================================================

        uv run python .\src\run_bnc.py ^
            --instance_file "%%~fF" ^
            --instance_type bobilib ^
            --lower_level general ^
            --projected 1 ^
            --cuts %%M ^
            --separation integer ^
            --time_lim %TIME_LIMIT% ^
            --max_cuts 20 ^
            --cplex_cuts -1 ^
            --only_root_node 0 ^
            --tolerance 1e-6 ^
            --verbose_level 0 ^
            --write_lps 0 ^
            --output_csv "!OUTPUT_CSV!"

        if errorlevel 1 (
            echo.
            echo ERROR: Run failed.
            echo Instance: %%~nF
            echo Method:   %%M
            echo Partial results remain in: !OUTPUT_CSV!
            popd
            exit /b 1
        )
    )
)

echo.
echo ============================================================
echo All N50 D50 experiments completed successfully.
echo Total runs: !RUN_COUNT!
echo Results:    !OUTPUT_CSV!
echo ============================================================

popd
endlocal
exit /b 0