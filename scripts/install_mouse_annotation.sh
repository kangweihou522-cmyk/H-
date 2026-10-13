#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
# Use Bioconductor 3.16 for the Debian/R 4.2 runtime; TLS and package MD5 remain enabled.
docker build --build-arg HTTP_PROXY --build-arg HTTPS_PROXY -f Dockerfile.annotation -t mouse-rnaseq:annotated .
docker tag mouse-rnaseq:annotated mouse-rnaseq:1
docker run --rm mouse-rnaseq:1 Rscript -e 'stopifnot(requireNamespace("org.Mm.eg.db",quietly=TRUE))'
