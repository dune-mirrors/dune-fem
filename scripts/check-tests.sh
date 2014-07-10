#!/bin/bash

if test $# -lt 1 ; then
  echo "Usage: $0 <dune-fem-dir> [MAKE_CHECK_FLAGS]"
  exit 1
fi

MAKE_CHECK_FLAGS=$2

WORKINGDIR=`pwd`
cd $1
FEMDIR=`pwd`

CHECKLOG=$WORKINGDIR/check-tests.out
TESTLOG=$WORKINGDIR/test.log
make -i check $MAKE_CHECK_FLAGS  &> $CHECKLOG

retvalue=0
warnings=`grep warning: $CHECKLOG | grep -v "default CommunicationManager is used" | grep -v "GRIDDIM not defined" | grep -v "No GRIDTYPE defined" | grep -v "Hdiv-Projection only working for polOrd = 1" | grep -v "YaspGrid does not provide a HierarchicIndexSet" | wc -l`
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

fails=$(grep "^FAIL:" $CHECKLOG | wc -l)
if test $fails -gt 0 ; then
  echo "Error: The following $fails tests have failed:"
  failedTests=$(grep "^FAIL:" $CHECKLOG | sed -r 's/FAIL: //g')
  echo $failedTests
  for file in $(grep "/test-suite.log" $CHECKLOG | sed -r 's/See //g'); do
    cat $file
  done &> $TESTLOG
  retvalue=1
fi

#if test x`grep "\\"All \\[\\[\\:digit\\:\\]\\]\\+ tests passed\\"" $CHECKLOG` == x ; then
#  exit 1
#fi

exit $retvalue
