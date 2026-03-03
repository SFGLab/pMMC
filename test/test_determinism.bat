@echo off
REM ============================================================
REM  test_determinism.bat - R-DET: Verify deterministic output
REM  Tests that fixed seed produces identical synthetic data
REM  (Uses 'generate' action which is fast and fully deterministic)
REM ============================================================
setlocal

set CUDAMMC=%~dp0..\build\Debug\pMMC.exe
if not exist "%CUDAMMC%" set CUDAMMC=%~dp0..\build\Release\pMMC.exe
if not exist "%CUDAMMC%" (
    echo [SKIP] pMMC.exe not found. Build first.
    exit /b 1
)

set OUTDIR=%TEMP%\cudammc_det_test
set ERRORS=0

echo ============================================================
echo  Determinism Test (R-DET)
echo ============================================================
echo.

REM Clean previous runs
rmdir /s /q "%OUTDIR%\run1" 2>nul
rmdir /s /q "%OUTDIR%\run2" 2>nul
mkdir "%OUTDIR%\run1" 2>nul
mkdir "%OUTDIR%\run2" 2>nul

echo [1/3] Generating synthetic data with seed 42 (run 1)...
"%CUDAMMC%" -a generate -o "%OUTDIR%/run1/" -j 42 -l 50 -m 1 -n det_test >"%OUTDIR%\run1\stdout.txt" 2>&1
if errorlevel 1 (
    echo [FAIL] Run 1 generate failed
    set /a ERRORS=ERRORS+1
    goto :results
)

echo [2/3] Generating synthetic data with seed 42 (run 2)...
"%CUDAMMC%" -a generate -o "%OUTDIR%/run2/" -j 42 -l 50 -m 1 -n det_test >"%OUTDIR%\run2\stdout.txt" 2>&1
if errorlevel 1 (
    echo [FAIL] Run 2 generate failed
    set /a ERRORS=ERRORS+1
    goto :results
)

REM Compare generated data files
echo.
echo [3/3] Comparing outputs...
set COMPARE_PASS=1

REM Compare anchors file
fc "%OUTDIR%\run1\anchors.txt" "%OUTDIR%\run2\anchors.txt" >nul 2>&1
if errorlevel 1 (
    echo [FAIL] anchors.txt differ between runs
    set COMPARE_PASS=0
    set /a ERRORS=ERRORS+1
) else (
    echo [PASS] anchors.txt identical
)

REM Compare clusters file
fc "%OUTDIR%\run1\clusters.txt" "%OUTDIR%\run2\clusters.txt" >nul 2>&1
if errorlevel 1 (
    echo [FAIL] clusters.txt differ between runs
    set COMPARE_PASS=0
    set /a ERRORS=ERRORS+1
) else (
    echo [PASS] clusters.txt identical
)

REM Compare singletons file
fc "%OUTDIR%\run1\singletons.txt" "%OUTDIR%\run2\singletons.txt" >nul 2>&1
if errorlevel 1 (
    echo [FAIL] singletons.txt differ between runs
    set COMPARE_PASS=0
    set /a ERRORS=ERRORS+1
) else (
    echo [PASS] singletons.txt identical
)

:results
echo.
echo ============================================================
if "%ERRORS%"=="0" (
    echo  DETERMINISM TEST PASSED
    exit /b 0
) else (
    echo  ERRORS: %ERRORS%
    exit /b 1
)
