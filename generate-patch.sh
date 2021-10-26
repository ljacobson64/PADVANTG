#!/bin/bash

base_dir=${PWD}
versions="3.0.3 3.2.0"

mkdir -pv patch
cd ${base_dir}/advantg-src

for version in ${versions}; do
  echo "Generating patch file for version ${version}"
  patch_file=${base_dir}/patch/${version}.patch
  diff -rN -U0 ${version}-vendor ${version}-perturb > ${patch_file}
  sed -i -e "s/.py\t.*/.py/" ${patch_file}
done
