# Build the Singularity/Apptainer image from chd-soybean.def (run on the production server).
# Needs internet access (Docker Hub, CRAN, Bioconductor, GitHub, PyPI); takes ~30-60 min,
# mostly compiling R packages. Ends with the image's %test (loads the R and Python packages).
# Then run the steps from the repo dir with the image, e.g.:
#   singularity exec --nv --bind /falafel $SIF bash step-3.sh   # --nv: GPU; --bind: data dirs
set -euo pipefail
cd "$(dirname "$0")" # run from the repo dir (chd-soybean.def copies requirements.txt from here)
source ./paths.sh # WORK

SIF=$WORK/chd-soybean.sif # the image (~10 GB or more)
# where the image is assembled during the build (~15-20 GB at its peak; emptied afterwards)
BUILD_TMPDIR=$WORK/apptainer-tmp
# downloaded base images (small; reused by rebuilds; empty with "singularity cache clean")
BUILD_CACHEDIR=$WORK/apptainer-cache
# --fakeroot builds without root if the admins enabled it; otherwise run this script with
# sudo and set BUILD_OPTS=""
BUILD_OPTS="--fakeroot"

# Apptainer or SingularityCE, whichever is installed (they read differently named variables)
if command -v apptainer > /dev/null; then
	SINGULARITY=apptainer
else
	SINGULARITY=singularity
fi
export APPTAINER_TMPDIR=$BUILD_TMPDIR SINGULARITY_TMPDIR=$BUILD_TMPDIR
export APPTAINER_CACHEDIR=$BUILD_CACHEDIR SINGULARITY_CACHEDIR=$BUILD_CACHEDIR
mkdir -p $BUILD_TMPDIR $BUILD_CACHEDIR

$SINGULARITY --version
$SINGULARITY build $BUILD_OPTS $SIF chd-soybean.def
rmdir $BUILD_TMPDIR 2> /dev/null || echo "note: $BUILD_TMPDIR is not empty (left over from a failed build?)"
echo "built $SIF"
