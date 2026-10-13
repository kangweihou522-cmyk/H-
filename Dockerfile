FROM debian:bookworm-slim
RUN apt-get update && apt-get install -y --no-install-recommends r-base-core r-bioc-deseq2 r-bioc-tximport r-bioc-go.db r-cran-jsonlite salmon && rm -rf /var/lib/apt/lists/*
WORKDIR /data
