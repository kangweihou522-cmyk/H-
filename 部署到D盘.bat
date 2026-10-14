@echo off
setlocal
powershell.exe -NoLogo -NoProfile -ExecutionPolicy Bypass -File "%~dp0windows\Workstation.ps1" -Action Install
set "RNASEQ_EXIT=%ERRORLEVEL%"
pause
exit /b %RNASEQ_EXIT%
