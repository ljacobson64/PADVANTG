#!/bin/bash

export PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/snap/bin
export PYTHONPATH=
export LIBRARY_PATH=
export LD_LIBRARY_PATH=
source ${HOME}/PADVANTG/advantg-install/3.2.0-perturb/advantg.rc
export MPIEXEC_PREFLAGS="--oversubscribe"
export MPIEXEC_POSTFLAGS="--oversubscribe"
#rm -rf model output fwd_solution adj_solution custom_output
advantg -v advantg.inp
