#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
rnaseq_venv="${RNASEQ_VENV:-$PWD/.venv}"
if [[ ! -x "$rnaseq_venv/bin/python" ]]; then echo '请先运行 bash install.sh'; exit 1; fi
printf '小鼠 RNA-seq 工作站：http://localhost:%s\n' "${PORT:-8765}"
exec "$rnaseq_venv/bin/python" "$PWD/app.py"
