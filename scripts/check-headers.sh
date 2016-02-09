#!/bin/bash

# check for parameter pointing to DUNE build directory
# ----------------------------------------------------

if test \( $# -lt 3 \) ; then
  echo "Usage: $0 <build-dir-base> <module-name> <MAKE_FLAGS>"
  exit 1
fi

# set up some variables
# ---------------------

WORKINGDIR=$(pwd)
cd $1
BUILDDIR=$(pwd)
MODULENAME=$2
MAKE_FLAGS=$3


errors=0

MODULES=""
for modctrl in $(ls -d */); do
  MODULES+=" $(basename $modctrl)"
done

for module in $MODULES; do
  cd $module
  echo
  echo "Checking headers in $module ..."
  CHECKLOG=$WORKINGDIR/headercheck-$module.log
  make headercheck -i $MAKE_FLAGS &> $CHECKLOG
  hc_errors=$(grep error: $CHECKLOG)
  if test -z "$hc_errors"; then
    rm $CHECKLOG
  else
    if test "x$module" = "x$MODULENAME" ; then
      echo "Error: headercheck for module $module failed (see headercheck-$module.log)"
      errors=$((errors+1))
    else
      echo "Warning: headercheck for module $module failed (see headercheck-$module.log)"
    fi
  fi
  cd $BUILDDIR
done


if test $errors -gt 0 ; then
  exit 1
else
  exit 0
fi
