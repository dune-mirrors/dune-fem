#!/bin/bash 

## start adaptation test with 4 MPI processes

# define executable 
export MPIEXEC_PROG=adaptation
# define executable 
export MPIEXEC_PARAMETER=parameter
# call mpi start script 
if [ -e ../../../scripts/mpiexec.sh ] ; then
  exec ../../../scripts/mpiexec.sh 4 $*
else
  if [ "$MPIEXEC" == "" ]; then
    mpiexec -np 4 ./adaptation
  else
    $MPIEXEC -np 4 ./adaptation
  fi
fi
