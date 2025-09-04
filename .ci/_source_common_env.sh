export CACHE_ROOT=$(pwd)/cache-ci
export PYTHONUSERBASE=$CACHE_ROOT/pyusrb
export CPATH=$SUNDIALS_ROOT/include:/usr/include/suitesparse:${CPATH:-}  # sunlinsol_klu.h includes "klu.h"
export CPLUS_INCLUDE_PATH=${Boost_ROOT}/include
export LIBRARY_PATH=$SUNDIALS_ROOT/lib
export LD_LIBRARY_PATH=$SUNDIALS_ROOT/lib

source $(compgen -G "/opt-3/cpython-v3.*-apt-deb/bin/activate")
