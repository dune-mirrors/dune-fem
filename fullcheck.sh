#!/bin/bash

errors=0
if test -e check.out ; then
  rm -f check.out
fi
touch check.out

if test -s check-dune-control.out ; then
  errors=$((errors+1))
  echo "**************************" >> check.out
  echo "Problem with dune-control:" >> check.out
  cat check-configure.out >> check.out
  echo "**************************" >> check.out
  exit 1
fi

./check-headers.sh \
| grep -v "Warning: Header is not included in Makefile.am: fem/io/visual/dx/dxdata.hh" \
| grep -v "Warning: Makefile.am not found in fem/operator/feop" \
| grep -v "Warning: Header is not included in Makefile.am: fem/quadrature/ugquadratures.hh" \
| grep -v "Warning: Makefile.am not found in fem/devel-suite/localdg/advectdiffusion/euler_mhd" \
>& check-headers.out
if test -s check-headers.out ; then
  errors=$((errors+1))
  echo "*********************" >> check.out
  echo "Problem with headers:" >> check.out
  cat check-headers.out >> check.out
  echo "*********************" >> check.out
fi

./check-doxygen.sh \
>& check-doxygen.out
if test -s check-doxygen.out ; then
  errors=$((errors+1))
  echo "***************************" >> check.out
  echo "Problem with documentation:" >> check.out
  cat check-doxygen.out >> check.out
  echo "***************************" >> check.out
fi
if [ $errors -gt 0 ] ; then
  echo "Number of errors: $errors"
  exit 1
fi
exit 0

