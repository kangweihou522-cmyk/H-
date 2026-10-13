@echo off
setlocal
if not defined WSL_DISTRO set WSL_DISTRO=Ubuntu-24.04
wsl -d %WSL_DISTRO% -- bash -lc "cd \"$(wslpath '%~dp0')\" && bash launch.sh"
pause
