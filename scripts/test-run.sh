#!/bin/bash

FEMDIR=$1
SCRITPSDIR="$FEMDIR/scripts"

errors=0
cd $SCRIPTSDIR
if test -e run.out ; then
  rm -f run.out
fi

if ! ./check-tests.sh $FEMDIR ; then
  errors=$((errors+1))
  echo "***********************" >> run.out
  echo "Problem with fem tests:" >> run.out
  cat check-tests.out >> run.out
  echo "***********************" >> run.out
else
  rm -f check-tests.out
fi

if [ $errors -gt 0 ] ; then
  echo "Number of errors: $errors"
  exit 1
fi



exit 0

