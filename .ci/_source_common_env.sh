export CACHE_ROOT=$(pwd)/cache-ci
export PYTHONUSERBASE=$CACHE_ROOT/pyusrb
export SUITESPARSE_ROOT=${SUITESPARSE_ROOT:-$(compgen -G "/opt-3/suitesparse-v*-release" | head -n1)}
# pykinsol still `#include`s "kinsol_spils.h", a header upstream SUNDIALS removed
# after the 6.x series, so it cannot build against the default $SUNDIALS_ROOT
# (7.8.0). It gets built against this older release instead, see run-01.
export PYKINSOL_SUNDIALS_ROOT=/opt-3/sundials-6.7.0-release
# pyodesys' native (JIT-compiled) backends #include "anyode/anyode.hpp", a
# header-only library pycvodes bundles as a build-time dependency but does not
# install; run-01 extracts it from pycvodes' own sdist into this directory.
export ANYODE_INCLUDE_ROOT=$CACHE_ROOT/anyode-include
export CPATH=$SUNDIALS_ROOT/include:$SUITESPARSE_ROOT/include/suitesparse:$ANYODE_INCLUDE_ROOT:${CPATH:-}  # sunlinsol_klu.h includes "klu.h"
export CPLUS_INCLUDE_PATH=${Boost_ROOT}/include
export LIBRARY_PATH=$SUNDIALS_ROOT/lib:$SUITESPARSE_ROOT/lib
export LD_LIBRARY_PATH=$SUNDIALS_ROOT/lib:$SUITESPARSE_ROOT/lib:$PYKINSOL_SUNDIALS_ROOT/lib
export PIP_FIND_LINKS=$CACHE_ROOT/wheelhouse  # pre-built pykinsol wheel, see run-01

source $(compgen -G "/opt-3/cpython-v3.*-apt-deb/bin/activate")
