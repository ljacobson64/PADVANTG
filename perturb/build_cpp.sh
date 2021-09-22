#!/bin/bash

set -e

jobs=`grep -c processor /proc/cpuinfo`
install_prefix=${PWD}
cnpy_dir=/opt/software_native/cnpy

cmake_string=
cmake_string+=" -DCNPY_DIR=${cnpy_dir}"
cmake_string+=" -DCMAKE_INSTALL_PREFIX=${install_prefix}"
cmake_string+=" -DCMAKE_BUILD_TYPE=Release"

rm -rf bld
mkdir bld
cd bld
cmake .. ${cmake_string}
make -j${jobs}
make install
cd ..
#rm -rf bld
