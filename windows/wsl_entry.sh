#!/usr/bin/env bash
set -euo pipefail
action="${1:?Missing action}"
root="${2:?Missing station root}"
program="${3:?Missing program path}"
rnaseq_venv="${4:-}"
port="${5:-8767}"
windows_root="${6:-}"

if [[ "$action" == prerequisites ]]; then
  if ! command -v python3 >/dev/null || ! python3 -c 'import venv, ensurepip' >/dev/null 2>&1 || ! command -v fc-list >/dev/null || [[ -z "$(fc-list :lang=zh 2>/dev/null)" ]]; then
    echo '首次准备 Ubuntu 的 Python 环境；sudo 如有提示，请输入你本机 Ubuntu 密码。'
    sudo apt-get update
    sudo apt-get install -y python3 python3-venv fonts-noto-cjk
  fi
  exit 0
fi

export RNASEQ_HOME="$root"
export RNASEQ_WINDOWS_HOME="$windows_root"
export RNASEQ_RUNS="$root/分析结果"
export RNASEQ_VENV="$rnaseq_venv"
export PORT="$port"
rnaseq_runtime="$(dirname "$rnaseq_venv")"
export MPLCONFIGDIR="$rnaseq_runtime/matplotlib"
export PIP_CACHE_DIR="$rnaseq_runtime/pip-cache"
export BUILDX_CONFIG="$rnaseq_runtime/buildx"
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-$rnaseq_runtime/cache}"
mkdir -p "$MPLCONFIGDIR" "$PIP_CACHE_DIR" "$BUILDX_CONFIG" "$XDG_CACHE_HOME" "$root/日志"
cd "$program"

case "$action" in
  install)
    bash install.sh
    ;;
  check)
    test -x "$rnaseq_venv/bin/python" || { echo '尚未安装 Python 环境，请运行“部署到D盘.bat”。'; exit 1; }
    command -v docker >/dev/null || { echo 'Ubuntu 无法找到 Docker；请启用 Docker Desktop 的 WSL Integration。'; exit 1; }
    docker info --format '{{.ServerVersion}}'
    docker image inspect mouse-rnaseq:2 >/dev/null
    "$rnaseq_venv/bin/python" -c 'import flask, pandas, numpy, scipy, matplotlib, adjustText, multiqc; print("Python dependencies: ready")'
    printf '程序：%s\n结果：%s\n网页：http://localhost:%s\n' "$program" "$RNASEQ_RUNS" "$port"
    ;;
  start)
    test -x "$rnaseq_venv/bin/python" || { echo '尚未安装依赖，请先运行“部署到D盘.bat”。'; exit 1; }
    docker image inspect mouse-rnaseq:2 >/dev/null || { echo '分析镜像不可用，请检查 Docker 或重新运行部署。'; exit 1; }
    exec 9>"$root/配置/server.lock"
    if ! flock -n 9; then echo '工作站已在启动或运行中。'; exit 0; fi
    state="$(python3 windows/station_runtime.py probe "$root" "$port")"
    mode="$(python3 -c 'import sys,json;print(json.load(sys.stdin)["state"])' <<< "$state")"
    if [[ "$mode" == running ]]; then echo '工作站已在运行。'; exit 0; fi
    if [[ "$mode" != offline ]]; then echo '端口被其他程序占用；未停止其他程序，请更换端口。'; exit 1; fi
    printf '转录组分析工作站\n程序：%s\n分析结果：%s\n网页：http://localhost:%s\n请保持此终端运行；按 Ctrl+C 关闭。\n' "$windows_root" "$RNASEQ_RUNS" "$port"
    exec bash launch.sh > >(tee -a "$root/日志/server.log") 2> >(tee -a "$root/日志/server.log" >&2)
    ;;
  *) echo "Unknown action: $action"; exit 2 ;;
esac
