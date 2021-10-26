#!/bin/bash

tally_list=$(cat tally_list.txt)
mkdir -p hdf5
cd hdf5
ln -snf ../forward/custom_output/xs.h5 advantg_xs.h5
ln -snf ../forward/fwd_solution/denovo-forward.inp.h5 advantg_fwd_inp.h5
ln -snf ../forward/fwd_solution/denovo-forward.out.h5 advantg_fwd_out.h5
for tally in ${tally_list}; do
  ln -snf ../adjoint-${tally}/adj_solution/denovo-adjoint.inp.h5 advantg_adj_${tally}_inp.h5
  ln -snf ../adjoint-${tally}/adj_solution/denovo-adjoint.out.h5 advantg_adj_${tally}_out.h5
done
