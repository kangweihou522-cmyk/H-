#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
command -v python3 >/dev/null || { echo '请先安装 python3 python3-venv'; exit 1; }
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
command -v docker >/dev/null || { echo '请先安装 Docker 并启动服务'; exit 1; }
docker info >/dev/null
build_args=(--build-arg HTTP_PROXY --build-arg HTTPS_PROXY)
# Managed cloud Docker DNS may not resolve the host's proxy alias. Reuse its resolved address.
proxy_name=$(python3 -c 'import os,urllib.parse;print(urllib.parse.urlparse(os.environ.get("HTTPS_PROXY",os.environ.get("HTTP_PROXY",""))).hostname or "")')
if [[ -n "$proxy_name" ]] && command -v getent >/dev/null; then
  proxy_addr=$(getent ahostsv4 "$proxy_name" | head -n 1 | awk '{print $1}') || true
  if [[ -n "${proxy_addr:-}" ]]; then build_args+=(--add-host "$proxy_name:$proxy_addr"); fi
fi
docker build "${build_args[@]}" -t mouse-rnaseq:core .
docker run --rm mouse-rnaseq:core Rscript -e 'stopifnot(all(sapply(c("DESeq2","tximport","GO.db","jsonlite"),requireNamespace,quietly=TRUE)))'
if docker image inspect mouse-rnaseq:annotated >/dev/null 2>&1; then
  docker tag mouse-rnaseq:annotated mouse-rnaseq:1
else
  docker tag mouse-rnaseq:core mouse-rnaseq:1
fi
echo '安装完成：bash launch.sh。小鼠注释库可用 scripts/install_mouse_annotation.sh 安装，或提供离线注释 TSV。'
