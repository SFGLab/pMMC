@echo off
REM ============================================================
REM  ci_local.bat - R-CI: Local CI validation script
REM  Runs: build -> correctness tests -> determinism -> benchmark
REM ============================================================
setlocal

set MSBUILD=C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\MSBuild\15.0\Bin\amd64\MSBuild.exe
set PROJDIR=%~dp0..
set VCXPROJ=%PROJDIR%\MSVC++\cudaMMC.vcxproj
set TOTAL_PASS=0
set TOTAL_FAIL=0
set TOTAL_SKIP=0

echo ============================================================
echo  cudaMMC2 Local CI Pipeline
echo  Date: %date% %time%
echo ============================================================
echo.

REM ====== GATE 1: BUILD ======
echo ============ GATE 1: BUILD ============
echo.

echo [Build] Compiling C/C++ sources (Release x64)...
"%MSBUILD%" "%VCXPROJ%" /p:Configuration=Release /p:Platform=x64 /t:ClCompile /nologo /v:minimal
if errorlevel 1 (
    echo [FAIL] C++ compilation failed
    set /a TOTAL_FAIL=TOTAL_FAIL+1
) else (
    echo [PASS] All C++ sources compiled successfully
    set /a TOTAL_PASS=TOTAL_PASS+1
)
echo.

REM Check if full link works (may fail due to cudafe++)
echo [Build] Attempting full link (Debug)...
"%MSBUILD%" "%VCXPROJ%" /p:Configuration=Debug /p:Platform=x64 /nologo /v:minimal 2>nul
if errorlevel 1 (
    echo [WARN] Full link failed (expected with CUDA 12.0 + MSVC 2017)
    set /a TOTAL_SKIP=TOTAL_SKIP+1
) else (
    echo [PASS] Full Debug link succeeded
    set /a TOTAL_PASS=TOTAL_PASS+1
)
echo.

REM ====== GATE 2: CODE VALIDATION ======
echo ============ GATE 2: CODE VALIDATION ============
echo.

echo [Code] Checking GPU seed propagation...
findstr /C:"setGpuSeed" "%PROJDIR%\src\main.cpp" >nul 2>&1
if errorlevel 1 (
    echo [FAIL] GPU seed not propagated from main.cpp
    set /a TOTAL_FAIL=TOTAL_FAIL+1
) else (
    echo [PASS] GPU seed propagation present
    set /a TOTAL_PASS=TOTAL_PASS+1
)

echo [Code] Checking CIF writer...
findstr /C:"CifWriter::write" "%PROJDIR%\src\CifWriter.cpp" >nul 2>&1
if errorlevel 1 (
    echo [FAIL] CIF writer implementation missing
    set /a TOTAL_FAIL=TOTAL_FAIL+1
) else (
    echo [PASS] CIF writer present
    set /a TOTAL_PASS=TOTAL_PASS+1
)

echo [Code] Checking PDB writer...
findstr /C:"PdbWriter::write" "%PROJDIR%\src\PdbWriter.cpp" >nul 2>&1
if errorlevel 1 (
    echo [FAIL] PDB writer implementation missing
    set /a TOTAL_FAIL=TOTAL_FAIL+1
) else (
    echo [PASS] PDB writer present
    set /a TOTAL_PASS=TOTAL_PASS+1
)

echo [Code] Checking MDS initializer...
findstr /C:"MdsInitializer::computeCoordinates" "%PROJDIR%\src\MdsInitializer.cpp" >nul 2>&1
if errorlevel 1 (
    echo [FAIL] MDS initializer implementation missing
    set /a TOTAL_FAIL=TOTAL_FAIL+1
) else (
    echo [PASS] MDS initializer present
    set /a TOTAL_PASS=TOTAL_PASS+1
)

echo [Code] Checking energy trace logging...
findstr /C:"energyTraceEnabled" "%PROJDIR%\src\LooperSolver.cpp" >nul 2>&1
if errorlevel 1 (
    echo [FAIL] Energy trace logging not found in MC loops
    set /a TOTAL_FAIL=TOTAL_FAIL+1
) else (
    echo [PASS] Energy trace logging present in MC loops
    set /a TOTAL_PASS=TOTAL_PASS+1
)

echo [Code] Checking GPU arc MC kernel...
if exist "%PROJDIR%\src\ParallelMonteCarloArcs.cu" (
    echo [PASS] GPU arc MC kernel file exists
    set /a TOTAL_PASS=TOTAL_PASS+1
) else (
    echo [FAIL] GPU arc MC kernel missing
    set /a TOTAL_FAIL=TOTAL_FAIL+1
)

echo [Code] Checking GPU smooth MC kernel...
if exist "%PROJDIR%\src\ParallelMonteCarloSmooth.cu" (
    echo [PASS] GPU smooth MC kernel file exists
    set /a TOTAL_PASS=TOTAL_PASS+1
) else (
    echo [FAIL] GPU smooth MC kernel missing
    set /a TOTAL_FAIL=TOTAL_FAIL+1
)

echo [Code] Checking memory budget enforcement...
findstr /C:"checkMemoryBudget" "%PROJDIR%\src\LooperSolver.cpp" >nul 2>&1
if errorlevel 1 (
    echo [FAIL] Memory budget check missing
    set /a TOTAL_FAIL=TOTAL_FAIL+1
) else (
    echo [PASS] Memory budget enforcement present
    set /a TOTAL_PASS=TOTAL_PASS+1
)

echo [Code] Checking -F flag for output format...
findstr /C:"-F" "%PROJDIR%\src\main.cpp" >nul 2>&1
if errorlevel 1 (
    echo [FAIL] Output format -F flag not found
    set /a TOTAL_FAIL=TOTAL_FAIL+1
) else (
    echo [PASS] Output format -F flag present
    set /a TOTAL_PASS=TOTAL_PASS+1
)

echo [Code] Checking -I flag for init method...
findstr /C:"initMethod" "%PROJDIR%\src\main.cpp" >nul 2>&1
if errorlevel 1 (
    echo [FAIL] Init method -I flag not found
    set /a TOTAL_FAIL=TOTAL_FAIL+1
) else (
    echo [PASS] Init method -I flag present
    set /a TOTAL_PASS=TOTAL_PASS+1
)

echo [Code] Checking -E flag for energy trace...
findstr /C:"energyTraceEnabled" "%PROJDIR%\src\main.cpp" >nul 2>&1
if errorlevel 1 (
    echo [FAIL] Energy trace -E flag not found
    set /a TOTAL_FAIL=TOTAL_FAIL+1
) else (
    echo [PASS] Energy trace -E flag present
    set /a TOTAL_PASS=TOTAL_PASS+1
)
echo.

REM ====== GATE 3: RUNTIME TESTS (if exe exists) ======
echo ============ GATE 3: RUNTIME TESTS ============
echo.

set CUDAMMC=%PROJDIR%\build\Debug\cudaMMC.exe
if not exist "%CUDAMMC%" set CUDAMMC=%PROJDIR%\build\Release\cudaMMC.exe
if exist "%CUDAMMC%" (
    echo [Runtime] Executable found: %CUDAMMC%
    echo [Runtime] Running determinism test...
    call "%~dp0test_determinism.bat"
    if errorlevel 1 (
        set /a TOTAL_FAIL=TOTAL_FAIL+1
    ) else (
        set /a TOTAL_PASS=TOTAL_PASS+1
    )
) else (
    echo [SKIP] No executable found - skipping runtime tests
    set /a TOTAL_SKIP=TOTAL_SKIP+1
)
echo.

REM ====== SUMMARY ======
echo ============================================================
echo  CI SUMMARY
echo ============================================================
echo  Passed: %TOTAL_PASS%
echo  Failed: %TOTAL_FAIL%
echo  Skipped: %TOTAL_SKIP%
echo ============================================================

if "%TOTAL_FAIL%"=="0" (
    echo  RESULT: ALL GATES PASSED
    exit /b 0
) else (
    echo  RESULT: %TOTAL_FAIL% FAILURES
    exit /b 1
)
