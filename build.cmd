@echo off
rem ===================================================================
rem  One-click build (MinGW g++)
rem  Eigen is header-only: tell the compiler where the folder that
rem  contains Eigen\Dense is with -I. This script finds it under include\.
rem ===================================================================
setlocal
cd /d "%~dp0"

set "EIGEN="
for /d %%D in ("include\eigen-*") do (
    if exist "%%D\Eigen\Dense" set "EIGEN=%%D"
)
rem also support Eigen\ placed directly under include\
if not defined EIGEN (
    if exist "include\Eigen\Dense" set "EIGEN=include"
)
if not defined EIGEN (
    echo [ERROR] Eigen headers not found. Expected: include\eigen-x.y.z\Eigen\Dense
    exit /b 1
)

echo [INFO] Eigen: %EIGEN%
echo [INFO] Compiling main.cpp -^> main.exe ...

rem -g -O0 : keep debug info so F5 debugging is easy to follow
g++ -std=c++17 -g -O0 -Wall -Wextra -I "%EIGEN%" main.cpp -o main.exe
if errorlevel 1 (
    echo.
    echo [FAILED] Compilation error, see messages above.
    exit /b 1
)

echo [OK] main.exe created
exit /b 0