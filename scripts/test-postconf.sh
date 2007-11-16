#!/bin/bash

FEMDIR=$1
SCRITPSDIR="$FEMDIR/scripts"
errors=0

cd $SCRITPSDIR
if test -e postconf.out ; then
  rm -f postconf.out
fi
touch postconf.out

./check-doxygen.sh $FEMDIR \
>& check-doxygen.out
if test -s check-doxygen.out ; then
  errors=$((errors+1))
  echo "***************************" >> postconf.out
  echo "Problem with documentation:" >> postconf.out
  cat check-doxygen.out >> postconf.out
  echo "***************************" >> postconf.out
else
  rm -r check-doxygen.out
fi

exit 0

