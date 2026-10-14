#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
command -v python3 >/dev/null || { echo '请先安装 python3 python3-venv'; exit 1; }
sha256sum -c vendor/SHA256SUMS
rnaseq_venv="${RNASEQ_VENV:-$PWD/.venv}"
python3 -m venv "$rnaseq_venv"
"$rnaseq_venv/bin/pip" install -r requirements.txt
command -v docker >/dev/null || { echo '请先安装 Docker 并启动服务'; exit 1; }
docker info >/dev/null
source scripts/build_options.sh
docker build "${build_args[@]}" -t mouse-rnaseq:2-core .
docker run --rm mouse-rnaseq:2-core Rscript -e 'stopifnot(all(sapply(c("DESeq2","tximport","GO.db","jsonlite","fgsea","GSVA","limma","WGCNA","decoupleR"),requireNamespace,quietly=TRUE)))'
if docker image inspect mouse-rnaseq:2-annotated >/dev/null 2>&1; then
  docker tag mouse-rnaseq:2-annotated mouse-rnaseq:2
else
  docker tag mouse-rnaseq:2-core mouse-rnaseq:2
fi
echo 'v2 安装完成：bash launch.sh。真实小鼠GO可用 scripts/install_mouse_annotation.sh 安装，或提供离线注释 TSV。'
