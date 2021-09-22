#!/bin/bash

# Quit if any errors occur
set -e

# Package versions
hdf5_version=1.10.7
silo_version=4.10.2
openmpi_version=4.1.1
lava_version=1.1.1
exnihilo_version=5.4.0
advantg_version=3.0.3

# Important paths
dist_dir=${HOME}/dist
dependency_dir=/opt/software_native
mcnp_exe=/opt/software_misc/MCNP/bin/mcnp5
DATAPATH=/opt/software_misc/MCNP/MCNP_DATA

# Miscellaneous environment variables
num_cpus=`grep -c processor /proc/cpuinfo`

# Empty the path environment variables
export LD_LIBRARY_PATH=
export LIBRARY_PATH=
export PYTHONPATH=

# Locations of dependencies
openmpi_dir=${dependency_dir}/openmpi-${openmpi_version}
hdf5_dir=${dependency_dir}/hdf5-${hdf5_version}
silo_dir=${dependency_dir}/silo-${silo_version}
lava_dir=${dependency_dir}/lava-${lava_version}

# Paths to MPI compilers
CC=${openmpi_dir}/bin/mpicc
CXX=${openmpi_dir}/bin/mpic++
FC=${openmpi_dir}/bin/mpifort

# Get ADVANTG source code
build_prefix=${PWD}/ADVANTG-${advantg_version}
install_prefix=${PWD}/ADVANTG-${advantg_version}
echo "Removing existing files"
rm -rf ${build_prefix}
mkdir -pv ${build_prefix}/bld
cd ${build_prefix}
tarball=advantg-${advantg_version}.tar.gz
echo "Unpacking ${tarball}"
tar -xzf ${dist_dir}/advantg/${tarball}

# Patch source code
cd advantg
echo "Applying patch"
patch -p1 < ../../patch/${advantg_version}.patch

# Install multigroup cross sections
echo "Installing mgxs"
mkdir -p ${install_prefix}/mgxs
cd ${install_prefix}/mgxs
tar -xzf ${dist_dir}/advantg/mgxs.tar.gz

# Set CMake settings
export CMAKE_PREFIX_PATH=${openmpi_dir}:${hdf5_dir}:${silo_dir}:${lava_dir}
cmake_string=
cmake_string+=" -DADVANTG_DEBUG=ON"
cmake_string+=" -DCMAKE_BUILD_TYPE=Release"
cmake_string+=" -DMCNP_EXECUTABLE=${mcnp_exe}"
cmake_string+=" -DSCALE_DATA_DIR="
cmake_string+=" -DANISNLIB_SEARCH_PATH=${install_prefix}/mgxs"
cmake_string+=" -DCMAKE_INSTALL_PREFIX=${install_prefix}"
cmake_string+=" -DUSE_OPENMP=OFF"
cmake_string+=" -DDENOVO_IS_PARALLEL=ON"
cmake_string+=" -DCMAKE_C_COMPILER=${CC}"
cmake_string+=" -DCMAKE_CXX_COMPILER=${CXX}"
cmake_string+=" -DCMAKE_Fortran_COMPILER=${FC}"

# Build and install ADVANTG
echo "Building and installing ADVANTG"
cd ${build_prefix}/bld
cmake ../advantg ${cmake_string}
make -j${num_cpus}
make -j${num_cpus} install
