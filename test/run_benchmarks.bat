@echo off
REM ============================================================
REM  run_benchmarks.bat - R-G4: Benchmark harness
REM  Reports per-phase timing and GPU vs CPU speedup
REM ============================================================
setlocal

set CUDAMMC=%~dp0..\build\Debug\pMMC.exe
if not exist "%CUDAMMC%" set CUDAMMC=%~dp0..\build\Release\pMMC.exe
if not exist "%CUDAMMC%" (
    echo [SKIP] pMMC.exe not found. Build first.
    exit /b 1
)

set OUTDIR=%TEMP%\cudammc_benchmark
echo ============================================================
echo  cudaMMC2 Benchmark Harness (R-G4)
echo ============================================================
echo.

mkdir "%OUTDIR%" 2>nul

REM Run the built-in benchmark action
echo Running benchmark (generates synthetic data, reconstructs, compares)...
echo.

echo [1/1] Benchmark with ensemble=5, seed=42...
"%CUDAMMC%" -a benchmark -o "%OUTDIR%/" -j 42 -m 5
echo.

echo ============================================================
echo  Benchmark Results
echo ============================================================
echo.

if exist "%OUTDIR%\benchmark_results.csv" (
    echo Results CSV:
    type "%OUTDIR%\benchmark_results.csv"
    echo.
) else (
    echo [INFO] No CSV results found. Check benchmark output above.
)

echo.
echo Output directory: %OUTDIR%
echo ============================================================
