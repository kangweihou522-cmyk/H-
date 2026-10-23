#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "$0")"
command -v python3 >/dev/null || { echo 'Ubuntu 缺少 python3，请先安装 python3 python3-venv。'; exit 1; }
exec python3 windows/ubuntu_deploy.py "$@"
