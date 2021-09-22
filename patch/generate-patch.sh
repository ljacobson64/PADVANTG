#!/bin/bash

versions="3.0.3"

echo "Generating patch file"
cd ..
mkdir -p patch
for version in ${versions}; do
  patch_file=patch/${version}.patch
  diff -rN -U0 src-${version}-vendor src-${version}-perturb > ${patch_file}
  sed -i -e "s/.cmake\t.*/.cmake/" ${patch_file}
  sed -i -e "s/.txt\t.*/.txt/"     ${patch_file}
  sed -i -e "s/.py\t.*/.py/"       ${patch_file}
done
