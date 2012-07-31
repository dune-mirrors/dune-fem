#!/bin/bash 

# This script is developed by Robert Kloefkorn (with some modifications by me)

# where should the program run
EXECHOST=localhost

# current path
PWD=`pwd`

# executable, change to your needs 
PROG=$PWD/l2projection_petscmatrix
# your specific parameters
#PREPARAMETER="b cginverseoperator.hh:278"
PREPARAMETER=
PARAMETER=

#env | grep OMPI
# rank and size info (OpenMPI)
RANK=$OMPI_COMM_WORLD_RANK
SIZE=$OMPI_COMM_WORLD_SIZE


if [ "$1" == "" ] ; then 
  if [ "$GDB" == "gdb" ] ; then 
    exec xterm -T "P[ $RANK ] - $PROG" -e gdb --command=run.gdb --args $PROG
  else 
    if [ "$FILES" == "outfiles" ] ; then 
      exec $PROG $PARAMETER 2> $PWD/mpierr.$RANK 1> $PWD/mpiout.$RANK
    else 
      exec $PROG $PARAMETER
    fi  
  fi
else
  if [ "$2" == "outfiles" ] ; then 
    export FILES="outfiles"
  fi  
  if [ "$2" == "gdb" ] ; then 
    export GDB="gdb"
    echo "Start in Debug Mode"
    echo "$PREPARAMETER" > run.gdb
    echo "run $PARAMETER" >> run.gdb 
    # start mpi program (will go into first if 
    nice mpiexec -n $1 -host $EXECHOST "$0" 
  else    
    # start mpi program (will go into first if 
    nice mpiexec -n $1 "$0"
  fi  
fi  
