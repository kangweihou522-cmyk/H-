@echo off
setlocal
powershell.exe -NoLogo -NoProfile -ExecutionPolicy Bypass -File "%~dp0windows\Workstation.ps1" -Action Start
set "RNASEQ_EXIT=%ERRORLEVEL%"
if not "%RNASEQ_EXIT%"=="0" pause
exit /b %RNASEQ_EXIT%
