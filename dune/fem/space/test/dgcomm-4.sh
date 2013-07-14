#!/bin/bash 

## start dgcomm test with 4 MPI processes

# define executable 
export MPIEXEC_PROG=dgcomm
# define executable 
export MPIEXEC_PARAMETER=parameter
# call mpi start script 
if [ -e ../../../../scripts/mpiexec.sh ] ; then
  exec ../../../../scripts/mpiexec.sh 4
else
  mpiexec -np 4 ./dgcomm
fi
