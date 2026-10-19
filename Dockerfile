# syntax=docker/dockerfile:1
# ChromatinHD soybean multiome pipeline: R environment (steps 0, 1a) and
# Python environment (steps 1b-5).
# Mirrors the setup recorded in packages_R.txt / packages_Python.txt.
# Same steps as chd-soybean.def (Singularity/Apptainer; keep the two in sync), which is
# what the production server builds (build-image.sh). Docker build: docker build -t chd-soybean .

# ubuntu:26.04 as of 2026-10-19, pinned by digest so rebuilds start from the same base
FROM ubuntu:26.04@sha256:c098726dd9b9733f3e885f8de0ef375206a311b0a317eadc40d347b1115fb944

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

# fail if apt gives another R than the one tested (Bioconductor 3.22 needs R 4.5)
RUN R_VERSION=$(Rscript -e 'cat(as.character(getRversion()))') \
    && if [ "$R_VERSION" != "4.5.2" ]; then echo "expected R 4.5.2, got R $R_VERSION" >&2; exit 1; fi

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

# Python 3.12 (uv-managed) in /opt/venv with the exact pins of requirements.txt
# (compiled from requirements.in; torch 2.14.1, CUDA 13 build from PyPI, chromatinhd 0.4.3
# built from source with the gcc installed above).
COPY --from=ghcr.io/astral-sh/uv:0.9.26 /uv /usr/local/bin/uv
ENV UV_PYTHON_INSTALL_DIR=/opt/uv-python \
    VIRTUAL_ENV=/opt/venv \
    PATH=/opt/venv/bin:$PATH
# /opt/requirements.txt stays in the image as a record of the pins it was built from
COPY requirements.txt /opt/requirements.txt
RUN uv venv /opt/venv --python 3.12.12 \
    && uv pip install --no-cache -r /opt/requirements.txt

WORKDIR /workspace
CMD ["bash"]
