#!/bin/bash

versions="3.0.3"

dist_dir=${HOME}/dist
cd ..
for version in ${versions}; do
  echo "Removing existing files"
  rm -rf src-${version}-vendor src-${version}-perturb
  tarball=${dist_dir}/advantg/advantg-${version}.tar.gz
  echo "Unpacking ${tarball}"
  tar -xzf ${tarball}
  mv -v advantg src-${version}-vendor
  echo "Copying src-${version}-vendor to src-${version}-perturb"
  cp -rp src-${version}-vendor src-${version}-perturb 
  echo "Applying patch"
  cd src-${version}-perturb
  patch -p1 < ../patch/${version}.patch
  cd ..
done
