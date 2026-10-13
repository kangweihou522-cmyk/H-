#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
source scripts/build_options.sh
if [[ -n "${CODEX_PROXY_CERT:-}" ]]; then
  build_args+=(--secret id=proxy_ca,src=/etc/ssl/certs/ca-certificates.crt)
fi
docker build "${build_args[@]}" -f Dockerfile.annotation -t mouse-rnaseq:2-annotated .
docker tag mouse-rnaseq:2-annotated mouse-rnaseq:2
docker run --rm mouse-rnaseq:2 Rscript -e 'stopifnot(requireNamespace("org.Mm.eg.db",quietly=TRUE))'
