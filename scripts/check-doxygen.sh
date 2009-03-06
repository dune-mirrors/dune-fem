#!/bin/bash

if test $# -lt 1 ; then
  echo "Usage: $0 <dune-fem-dir>"
  exit 1
fi

cd $1
FEMDIR=`pwd`

errors=0

DOXYGEN_OUT="$FEMDIR/doc/doxygen/doxygen.out"
if test -s $DOXYGEN_OUT ; then
  errors=$((errors+1))
  echo "Warning: Doxygen output not empty (see $DOXYGEN_OUT)!"
fi

DOXYGEN_FORMULALOG="$FEMDIR/doc/doxygen/html/_formulas.log"
if test -s $DOXYGEN_FORMULALOG ; then
  errors=$((errors+1))
  echo "Warning: Problem in documentation with latex code (see $DOXYGEN_FORMULALOG)!"
fi

if test $errors -gt 0 ; then
  exit 1
else
  exit 0
fi
