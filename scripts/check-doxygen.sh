#!/bin/bash

if test $# -lt 1 ; then
  echo "Usage: $0 <dune-fem-dir>"
  exit 1
fi

cd $1

errors=0
if test -s ./doc/doxygen/doxygen.out ; then
  errors=$((errors+1))
  echo "Warning: Doxygen output not empty!"
fi
if test -s ./doc/doxygen/html/_formulas.log ; then
  errors=$((errors+1))
  echo "Warning: Problem in documentation with latex code - _formulars.log not empty!"
fi

if test $errors -gt 0 ; then
  exit 1
else
  exit 0
fi
