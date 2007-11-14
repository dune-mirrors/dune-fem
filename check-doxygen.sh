#!/bin/bash

if test -s ./doc/doxygen/doxygen.out ; then
  echo "Warning: Doxygen output not empty!"
fi
if test -s ./doc/doxygen/html/_formulas.log ; then
  echo "Warning: Problem in documentation with latex code - _formulars.log not empty!" 
fi

