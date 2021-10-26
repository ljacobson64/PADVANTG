#!/bin/bash

set -e

base_dir=${PWD}
versions="3.0.3 3.2.0"
src_types="vendor perturb"

rm    -rf advantg-src advantg-install
mkdir -pv advantg-src advantg-install

for version in ${versions}; do
  advantg_dir=/opt/software/ADVANTG-${version}

  # Location of relevant source files within install directory
  if [ ${version} == "3.0.3" ]; then
    subdirs="python"
  elif [ ${version} == "3.2.0" ]; then
    subdirs="packages advantg python"
  else
    echo "Invalid version: ${version}"
    exit 0
  fi
  src_dir=${advantg_dir}
  for subdir in ${subdirs}; do
    src_dir=${src_dir}/${subdir}
  done
  src_dir=${src_dir}/advantg

  # Set up source directories
  cd ${base_dir}/advantg-src
  cp -rp ${src_dir} ${version}-vendor
  cp -rp ${src_dir} ${version}-perturb
  cd ${version}-perturb
  patch -p1 < ${base_dir}/patch/${version}.patch

  # Set up install directories
  cd ${base_dir}/advantg-install
  mkdir -pv ${version}-vendor
  cd ${version}-vendor
  src_dir=${advantg_dir}
  for subdir in ${subdirs}; do
    for f in $(ls ${src_dir} -I ${subdir}); do
      ln -sv ${src_dir}/${f} .
    done
    src_dir=${src_dir}/${subdir}
    mkdir -pv ${subdir}
    cd ${subdir}
  done
  cd ${base_dir}/advantg-install
  cp -rp ${version}-vendor ${version}-perturb
  for src_type in ${src_types}; do
    cd ${base_dir}/advantg-install
    install_dir=${version}-${src_type}
    for subdir in ${subdirs}; do
      install_dir=${install_dir}/${subdir}
    done
    cd ${install_dir}
    ln -sv ${base_dir}/advantg-src/${version}-${src_type} advantg
  done

  # Fix advantg.rc
  if [ ${version} == "3.0.3" ]; then continue; fi
  for src_type in ${src_types}; do
    cd ${base_dir}/advantg-install/${version}-${src_type}
    rm -fv advantg.rc
    cp -pv ${advantg_dir}/advantg.rc .
    dir1=`echo ${advantg_dir}                                     | sed 's/\//\\\\\//g'`
    dir2=`echo ${base_dir}/advantg-install/${version}-${src_type} | sed 's/\//\\\\\//g'`
    sed -i "s/ADVANTG=${dir1}/ADVANTG=${dir2}/" advantg.rc
  done
done
