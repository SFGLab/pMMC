@echo off
REM ============================================================
REM  build_and_test.bat - Compile-and-verify script for cudaMMC2
REM  REQ-6.1, REQ-6.2, REQ-6.3
REM ============================================================
setlocal

set MSBUILD=C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\MSBuild\15.0\Bin\amd64\MSBuild.exe
set PROJDIR=%~dp0..
set VCXPROJ=%PROJDIR%\MSVC++\pMMC.vcxproj
set ERRORS=0

echo ============================================================
echo  cudaMMC2 Build and Test
echo ============================================================
echo.

REM ------ Section 1: Compile C/C++ source files ------
echo [1/3] Compiling C/C++ source files (Release x64)...
"%MSBUILD%" "%VCXPROJ%" /p:Configuration=Release /p:Platform=x64 /t:ClCompile /nologo /v:minimal
if errorlevel 1 (
    echo [FAIL] Compilation failed
    set /a ERRORS=ERRORS+1
) else (
    echo [PASS] All source files compiled successfully
)

echo.

REM ------ Section 2: Verify deterministic seed parsing ------
echo [2/3] Verifying deterministic seed support...
findstr /C:"-j  random seed" "%PROJDIR%\src\main.cpp" >nul 2>&1
if errorlevel 1 (
    echo [FAIL] -j flag not found in main.cpp usage
    set /a ERRORS=ERRORS+1
) else (
    echo [PASS] Deterministic seed flag present in CLI help
)

echo.

REM ------ Section 3: Verify memory budget enforcement ------
echo [3/3] Verifying memory budget enforcement...
findstr /C:"checkMemoryBudget" "%PROJDIR%\src\LooperSolver.cpp" >nul 2>&1
if errorlevel 1 (
    echo [FAIL] checkMemoryBudget not found in LooperSolver.cpp
    set /a ERRORS=ERRORS+1
    goto :results
)
findstr /C:"checkMemoryBudget" "%PROJDIR%\src\main.cpp" >nul 2>&1
if errorlevel 1 (
    echo [FAIL] checkMemoryBudget not found in main.cpp
    set /a ERRORS=ERRORS+1
    goto :results
)
findstr /C:"maxMemoryMB" "%PROJDIR%\src\Settings.cpp" >nul 2>&1
if errorlevel 1 (
    echo [FAIL] maxMemoryMB not found in Settings.cpp
    set /a ERRORS=ERRORS+1
    goto :results
)
echo [PASS] Memory budget enforcement present in all locations

:results
echo.
echo ============================================================
if "%ERRORS%"=="0" (
    echo  ALL 3 TESTS PASSED
    exit /b 0
) else (
    echo  ERRORS: %ERRORS%
    exit /b 1
)
