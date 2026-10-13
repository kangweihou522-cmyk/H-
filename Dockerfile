FROM debian:bookworm-slim@sha256:7c7b2c966bc9ee8cedfeef67e0e279108992c77681fa595db4a9d65c06ccc587
RUN apt-get update && apt-get install -y --no-install-recommends r-base-core r-base-dev r-bioc-deseq2 r-bioc-tximport r-bioc-go.db r-cran-jsonlite r-cran-data.table r-cran-cowplot r-cran-fastmatch r-bioc-gsva r-bioc-limma r-cran-wgcna r-cran-broom r-cran-tidyr r-cran-purrr salmon fastp && rm -rf /var/lib/apt/lists/*
COPY vendor/fgsea-1.24.0.tar.gz /tmp/fgsea.tar.gz
RUN mkdir /tmp/fgsea && tar -xzf /tmp/fgsea.tar.gz -C /tmp/fgsea --strip-components=1 && R CMD INSTALL /tmp/fgsea && rm -rf /tmp/fgsea /tmp/fgsea.tar.gz
COPY vendor/decoupleR-2.4.0.tar.gz /tmp/decoupleR.tar.gz
RUN mkdir /tmp/decoupleR && tar -xzf /tmp/decoupleR.tar.gz -C /tmp/decoupleR --strip-components=1 && R CMD INSTALL /tmp/decoupleR && rm -rf /tmp/decoupleR /tmp/decoupleR.tar.gz
WORKDIR /data
