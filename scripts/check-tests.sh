#!/bin/bash

if test $# -lt 1 ; then
  echo "Usage: $0 <dune-fem-dir>"
  exit 1
fi

WORKINGDIR=`pwd`
cd $1
FEMDIR=`pwd`

CHECKLOG=$WORKINGDIR/check-tests.out
make -i check &> $CHECKLOG

warnings=`grep warning: $CHECKLOG | wc -l`
if test $warnings -gt 0 ; then
  echo "Warning: $warnings compiler warnings occurred."
fi
errors=`grep error: $CHECKLOG | wc -l`
if test $errors -gt 0 ; then
  echo "Error: $errors compile time errors occurred."
  exit 1
fi

#if test x`grep "\\"All \\[\\[\\:digit\\:\\]\\]\\+ tests passed\\"" $CHECKLOG` == x ; then
#  exit 1
#fi

exit 0
