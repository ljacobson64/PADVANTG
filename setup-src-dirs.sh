#!/bin/bash

versions="3.0.3 3.2.0"
dist_dir=${HOME}/dist
for version in ${versions}; do
  echo "Removing existing files"
  rm -rf src-${version}-vendor src-${version}-perturb
  tarball=${dist_dir}/advantg/advantg-${version}.tar.gz
  echo "Unpacking ${tarball}"
  tar -xzf ${tarball}
  if [ "${version}" == "3.0.3" ]; then
    mv -v advantg src-${version}-vendor
  else
    mv -v advantg-${version} src-${version}-vendor
  fi
  echo "Copying src-${version}-vendor to src-${version}-perturb"
  cp -rp src-${version}-vendor src-${version}-perturb 
  echo "Applying patch"
  cd src-${version}-perturb
  patch -p1 < ../patch/${version}.patch
  cd ..
done
