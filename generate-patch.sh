#!/bin/bash

versions="3.0.3 3.2.0"
mkdir -pv patch
for version in ${versions}; do
  echo "Generating patch file for version ${version}"
  patch_file=patch/${version}.patch
  diff -rN -U0 src-${version}-vendor src-${version}-perturb > ${patch_file}
  sed -i -e "s/.cmake\t.*/.cmake/" ${patch_file}
  sed -i -e "s/.txt\t.*/.txt/"     ${patch_file}
  sed -i -e "s/.py\t.*/.py/"       ${patch_file}
done
