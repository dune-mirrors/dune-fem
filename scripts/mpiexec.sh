#!/bin/bash 

# where should the program run
EXECHOST=localhost

# current path
PWD=`pwd`

# your specific parameters, use export PARAMETER=...
if [ "$MPIEXEC_PARAMETER" == "" ]; then 
  MPIEXEC_PARAMETER="parameter"
fi

# executable to execute 
if [ "$MPIEXEC_PROG" == "" ]; then 
  echo "no executable given, export MPIEXEC_PROG='your program'"
  exit 1
fi  

# executable, change to your needs 
PROG=$PWD/$MPIEXEC_PROG 

#env | grep OMPI
# rank and size info (OpenMPI)
RANK=$OMPI_COMM_WORLD_RANK
SIZE=$OMPI_COMM_WORLD_SIZE

if [ "$1" == "" ] ; then 
  if [ "$GDB" == "gdb" ] ; then 
    exec xterm -T "$PROG: P[ $RANK ]" -e gdb --command=run.gdb --args $PROG
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
    echo "run $PARAMETER" > run.gdb 
    # start mpi program (will go into first if 
    nice mpiexec -n $1 -host $EXECHOST "$0" 
  else    
    # start mpi program (will go into first if 
    nice mpiexec -n $1 "$0"
  fi  
fi  
