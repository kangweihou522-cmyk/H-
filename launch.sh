#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
if [[ ! -x .venv/bin/python ]]; then echo '请先运行 bash install.sh'; exit 1; fi
printf '小鼠 RNA-seq 工作站：http://localhost:%s\n' "${PORT:-8765}"
exec .venv/bin/python app.py
