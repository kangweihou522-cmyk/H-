#!/usr/bin/env bash
set -euo pipefail
here="$(cd -- "$(dirname -- "$0")" && pwd)"
rnaseq_pwsh="${RNASEQ_PWSH:-pwsh}"
shim="$(mktemp -d /tmp/rnaseq-wsl-native.XXXXXXXX)"
trap 'rm -rf -- "$shim"' EXIT
cp "$here/fixtures/wsl_shim.py" "$shim/wsl.exe"
chmod +x "$shim/wsl.exe"
mkdir "$shim/second"
cp "$shim/wsl.exe" "$shim/second/wsl.exe"
printf ABC | PATH="$shim:$shim/second:$PATH" "$rnaseq_pwsh" -NoLogo -NoProfile -File "$here/check_windows_native.ps1"
