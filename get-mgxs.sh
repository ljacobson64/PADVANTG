#!/bin/bash

dist_dir=${HOME}/dist
mkdir -pv mgxs
cd mgxs
echo "Unpacking mgxs"
tar -xzf ${dist_dir}/advantg/mgxs.tar.gz
echo "Unpacking mgxs-more"
tar -xzf ${dist_dir}/advantg/mgxs-more.tar.gz
