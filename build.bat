@echo off
REM ============================================================
REM  build.bat — One-command build for pMMC (cudaMMC)
REM  Platform: Windows 10/11 with MSVC 2017 + CUDA 12.0
REM  Usage: cd pMMC\cudaMMC2 && build.bat [Release|Debug]
REM ============================================================

SET CONFIG=%1
IF "%CONFIG%"=="" SET CONFIG=Debug

SET MSBUILD="C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\MSBuild\15.0\Bin\amd64\MSBuild.exe"

IF NOT EXIST %MSBUILD% (
    echo ERROR: MSBuild not found at expected path.
    echo Please install Visual Studio 2017 with C++ and CUDA workloads.
    exit /b 1
)

echo Building pMMC (cudaMMC) in %CONFIG% configuration...
echo.

%MSBUILD% MSVC++\pMMC.vcxproj ^
    /p:Configuration=%CONFIG% ^
    /p:Platform=x64 ^
    /m ^
    /nologo ^
    /verbosity:minimal

IF %ERRORLEVEL% EQU 0 (
    echo.
    echo Build SUCCEEDED.
    echo Binary: build\%CONFIG%\pMMC.exe
) ELSE (
    echo.
    echo Build FAILED with error code %ERRORLEVEL%.
    echo NOTE: If Release build fails on ParallelMonteCarloHeatmap.cu,
    echo       this is a known cudafe++ issue with CUDA 12.0 + MSVC 2017.
    echo       Use Debug build: build.bat Debug
    exit /b %ERRORLEVEL%
)
