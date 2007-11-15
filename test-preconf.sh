#!/bin/bash

errors=0
if test -e preconf.out ; then
  rm -f preconf.out
fi
touch preconf.out

./check-headers.sh \
| grep -v "Warning: Header is not included in Makefile.am: fem/io/visual/dx/dxdata.hh" \
| grep -v "Warning: Makefile.am not found in fem/operator/feop" \
| grep -v "Warning: Header is not included in Makefile.am: fem/quadrature/ugquadratures.hh" \
| grep -v "Warning: Makefile.am not found in fem/devel-suite/localdg/advectdiffusion/euler_mhd" \
>& check-headers.out
if test -s check-headers.out ; then
  errors=$((errors+1))
  echo "*********************" >> preconf.out
  echo "Problem with headers:" >> preconf.out
  cat check-headers.out >> preconf.out
  echo "*********************" >> preconf.out
fi

exit 0

