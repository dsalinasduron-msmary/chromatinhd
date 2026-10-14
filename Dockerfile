# syntax=docker/dockerfile:1
# ChromatinHD soybean multiome pipeline - R environment (steps 0-1).
# Python environment (steps 2-5) to be added.
# Mirrors the setup recorded in packages_R.txt. For Singularity/Apptainer:
#   docker build -t chd-soybean .
#   singularity build chd-soybean.sif docker-daemon://chd-soybean:latest

FROM ubuntu:26.04

ENV DEBIAN_FRONTEND=noninteractive \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    TZ=Etc/UTC

# R 4.5.2 from Ubuntu, plus system libraries needed to build packages from source
RUN apt-get update && apt-get install -y --no-install-recommends \
        r-base r-base-dev \
        build-essential gfortran cmake git ca-certificates curl \
        libhdf5-dev \
        libcurl4-openssl-dev libssl-dev libxml2-dev libgit2-dev \
        zlib1g-dev libbz2-dev liblzma-dev \
        libglpk-dev libgsl-dev \
        libfontconfig1-dev libfreetype-dev libharfbuzz-dev libfribidi-dev \
        libpng-dev libtiff-dev libjpeg-dev \
    && rm -rf /var/lib/apt/lists/*

# CRAN frozen at a dated snapshot so rebuilds get the same versions
# (Seurat 5.6.0, Signac 1.17.1, ...); Bioconductor pinned to 3.22.
ARG CRAN_SNAPSHOT=https://packagemanager.posit.co/cran/2026-10-14
ARG BIOC_VERSION=3.22
ARG MUDATASEURAT_REF=e34e9082a3da029e6b21f7a93a193a2cc0d69245

RUN cat > /tmp/install_packages.R <<EOF
options(repos = c(CRAN = "${CRAN_SNAPSHOT}"), Ncpus = parallel::detectCores())
install.packages(c("remotes", "BiocManager"))
BiocManager::install(version = "${BIOC_VERSION}", update = FALSE, ask = FALSE)
# BiocManager resolves both CRAN and Bioconductor dependencies
# (Signac needs Bioconductor's sparseMatrixStats).
BiocManager::install(c(
    "dplyr", "readr", "magrittr", "stringr",
    "hdf5r", "Seurat", "SeuratObject", "Signac",
    "GenomicRanges", "GenomeInfoDb", "IRanges", "Rsamtools",
    "S4Vectors", "BiocGenerics", "sparseMatrixStats"
), update = FALSE, ask = FALSE)
# WriteH5MU comes from MuDataSeurat (GitHub only), not scCustomize (not needed)
remotes::install_github("PMBio/MuDataSeurat@${MUDATASEURAT_REF}", upgrade = "never")
pk <- c("dplyr", "readr", "magrittr", "stringr", "hdf5r", "Seurat",
        "SeuratObject", "Signac", "MuDataSeurat")
missing <- pk[!vapply(pk, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("failed to install: ", paste(missing, collapse = ", "))
EOF
RUN Rscript /tmp/install_packages.R && rm /tmp/install_packages.R

WORKDIR /workspace
CMD ["bash"]
