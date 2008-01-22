#!/bin/bash

if test x$SLURM_NPROCS == x ; then
  echo 'This script shall be used with slurm. For example use the following command:'
  echo 'srun -n <number of processes> --batch parallel.sh'
  exit 1
fi

export OMPI_MCA_mpi_yield_when_idle=1
mpiexec -np $SLURM_NPROCS nice -n 19 ./parallel 8
