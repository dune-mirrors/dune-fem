#!/bin/bash 

## start dgcomm test with 4 MPI processes

# define executable 
export MPIEXEC_PROG=dgcomm
# define executable 
export MPIEXEC_PARAMETER=parameter
# call mpi start script 
exec ../../../../scripts/mpiexec.sh 4
