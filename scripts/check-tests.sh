#!/bin/bash

FEMDIR=$1
SCRITPSDIR="$FEMDIR/scripts"

cd $FEMDIR
make check 2>&1 > $SCRIPTSDIR/check-tests.out
cd $SCRIPTSDIR
TEST=`grep "All [[:digit:]] tests passed" check-tests.out`
if test "x$TEST" == "x" ; then
  exit 1
fi

exit 0
