#!/usr/bin/env bash
# Source this helper from installers; registry credentials remain unchanged.
build_args=(--build-arg HTTP_PROXY --build-arg HTTPS_PROXY)
proxy_name=$(python3 -c 'import os,urllib.parse;print(urllib.parse.urlparse(os.environ.get("HTTPS_PROXY",os.environ.get("HTTP_PROXY",""))).hostname or "")')
if [[ -n "$proxy_name" ]] && command -v getent >/dev/null; then
  proxy_addr=$(getent ahostsv4 "$proxy_name" | head -n 1 | awk '{print $1}') || true
  if [[ -n "${proxy_addr:-}" ]]; then build_args+=(--add-host "$proxy_name:$proxy_addr"); fi
fi
