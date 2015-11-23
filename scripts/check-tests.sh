#!/bin/bash

if test $# -lt 1 ; then
  echo "Usage: $0 <dune-base-dir> [MAKE_CHECK_FLAGS]"
  exit 1
fi

MAKE_CHECK_FLAGS=$2

WORKINGDIR=`pwd`
cd $1
BUILDDIR=`pwd`

echo
echo "Performing make build_tests && make test in $(basename $BUILDDIR) ..."

CHECKLOG=$WORKINGDIR/check-tests.log
make -i build_tests $MAKE_CHECK_FLAGS  &> $CHECKLOG
make test &> /dev/null

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

# if a test has faild this file is writen
if test -e $BUILDDIR/Testing/Temporary/LastTestsFailed.log ; then
  retvalue=1
fi

# copy files with test logs and timings
mv $BUILDDIR/Testing/Temporary/* $WORKINGDIR/

exit $retvalue
