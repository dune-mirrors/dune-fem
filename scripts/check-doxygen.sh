#!/bin/bash

FEMDIR=$1
SCRITPSDIR="$FEMDIR/scripts"

errors=0
cd $FEMDIR
if test -s ./doc/doxygen/doxygen.out ; then
  errors=$((errors+1))
  echo "Warning: Doxygen output not empty!" >> $SCRIPTSDIR/check-doxygen.out
fi
if test -s ./doc/doxygen/html/_formulas.log ; then
  errors=$((errors+1))
  echo "Warning: Problem in documentation with latex code - _formulars.log not empty!" >> $SCRIPTSDIR/check-doxygen.out
fi

if [ $errors -gt 0 ] ; then
  echo "Number of errors: $errors"
  exit 1
fi

exit 0
