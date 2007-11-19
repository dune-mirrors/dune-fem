#!/bin/bash

if test $# -lt 1 ; then
  echo "Usage: $0 <dune-fem-dir>"
  exit 1
fi

WORKINGDIR=`pwd`
cd $1
FEMDIR=`pwd`

CHECKLOG=$WORKINGDIR/check-tests.out
make check &> $CHECKLOG
if test x`grep "All [[:digit:]]\+ tests passed" $CHECKLOG` == x ; then
  exit 1
else
  exit 0
fi
