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

retvalue=0
warnings=`grep warning: $CHECKLOG | wc -l`
if test $warnings -gt 0 ; then
  echo "Warning: $warnings compiler warnings occurred."
fi
errors=`grep error: $CHECKLOG | wc -l`
if test $errors -gt 0 ; then
  echo "Error: $errors compile time errors occurred."
  retvalue=1
fi
urefs=`grep ": undefined reference" $CHECKLOG | wc -l`
if test $urefs -gt 0 ; then
  echo "Error: $urefs undefined linker references occurred."
  retvalue=1
fi

#if test x`grep "\\"All \\[\\[\\:digit\\:\\]\\]\\+ tests passed\\"" $CHECKLOG` == x ; then
#  exit 1
#fi

exit $retvalue
