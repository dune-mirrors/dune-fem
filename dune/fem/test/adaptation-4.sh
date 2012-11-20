#!/bin/bash 

## start adaptation test with 4 MPI processes

# define executable 
export MPIEXEC_PROG=adaptation
# define executable 
export MPIEXEC_PARAMETER=parameter
# call mpi start script 
exec ../../../scripts/mpiexec.sh 4
