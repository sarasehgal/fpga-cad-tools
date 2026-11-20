#!/bin/bash
# builds just the parts of suitesparse we need (umfpack + amd + config) as static libs, no blas/cholmod needed.
# puts them in ./suitesparse/{include,lib}. only hv to run it once:  bash build_umfpack.sh [path to SuiteSparse-5.10.1.tar.gz]
cd "$(dirname "$0")"
TAR=${1:-../solver/SuiteSparse-5.10.1.tar.gz}
[ -f "$TAR" ] || TAR=solver/SuiteSparse-5.10.1.tar.gz
[ -f "$TAR" ] || TAR=SuiteSparse-5.10.1.tar.gz
[ -f "$TAR" ] || { echo "cant find SuiteSparse-5.10.1.tar.gz, pass the path to it"; exit 1; }
set -e
rm -rf ss_build && mkdir ss_build
tar xzf "$TAR" -C ss_build
SRC=ss_build/SuiteSparse-5.10.1
CFG="UMFPACK_CONFIG=-DNCHOLMOD -DNBLAS"
make -C $SRC/SuiteSparse_config static CFOPENMP= -j4 >/dev/null
make -C $SRC/AMD/Lib static CFOPENMP= -j4 >/dev/null
make -C $SRC/UMFPACK/Lib static CFOPENMP= "$CFG" -j4 >/dev/null
mkdir -p suitesparse/include suitesparse/lib
cp $SRC/SuiteSparse_config/SuiteSparse_config.h $SRC/AMD/Include/amd.h $SRC/UMFPACK/Include/*.h suitesparse/include/
cp $SRC/SuiteSparse_config/libsuitesparseconfig.a $SRC/AMD/Lib/libamd.a $SRC/UMFPACK/Lib/libumfpack.a suitesparse/lib/
rm -rf ss_build
echo "done -> suitesparse/lib: $(ls suitesparse/lib)"
