# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.50)
AC_INIT(dune-fem, 0.0, dune@mathematik.uni-freiburg.de)
AM_INIT_AUTOMAKE(dune-fem, 0.0, dune@mathematik.uni-freiburg.de)

# check if we're in the right directory
AC_CONFIG_SRCDIR([fem/space/dgspace/dgspace.hh])
# which header-file is created
AM_CONFIG_HEADER(config.h)

# we need no more than the standard DUNE-stuff
DUNE_CHECK_ALL_M([dune-common])
DUNE_CHECK_ALL_M([dune-grid])

DUNE_SUMMARY_ALL

if test x$HAVE_DUNE_COMMON != x1 ; then
  AC_MSG_ERROR([Can't work without the dune-common module. Maybe you have to supply your dune-common directory as --with-dunecommon=dir])
fi

if test x$HAVE_DUNE_GRID != x1 ; then
  AC_MSG_ERROR([Can't work without the dune-grid module. Maybe you have to supply your dune-common directory as --with-dunegrid=dir])
fi

# preset variable to path such that #include <dune/...> works
AC_SUBST(AM_CPPFLAGS, '$(DUNE_CPPFLAGS) -I$(top_srcdir)')
AC_SUBST(AM_LDFLAGS, '$(DUNE_LDFLAGS) $(DUNE_LIBS)')

# implicitly set the Dune-flags everywhere
CPPFLAGS="$DUNE_CPPFLAGS"
LDFLAGS="$DUNE_LDFLAGS"
LIBS="$DUNE_LIBS"

# generate Makefiles
AC_CONFIG_FILES([Makefile
     fem/Makefile
     fem/examples/Makefile
     fem/examples/poisson/Makefile
     fem/examples/elliptic/Makefile
     fem/examples/transport/Makefile
     fem/examples/vectorial/Makefile
     fem/examples/l2projection/Makefile
     fem/examples/localdg/Makefile
     fem/examples/localdg/advectdiffusion/Makefile
     fem/examples/pass/Makefile
     fem/pass/Makefile
     fem/pass/test/Makefile
     fem/pass/test/dgellipt/Makefile
     fem/space/Makefile
     fem/space/test/Makefile
     fem/space/basefunctions/Makefile
     fem/space/dgspace/Makefile
     fem/space/lagrangespace/Makefile
     fem/space/combinedspace/Makefile
     fem/space/common/Makefile
     fem/quadrature/Makefile
     fem/quadrature/test/Makefile
     fem/quadrature/caching/Makefile
     fem/quadrature/caching/test/Makefile
     fem/discretefunction/Makefile
     fem/discretefunction/common/Makefile
     fem/discretefunction/adaptivefunction/Makefile
     fem/discretefunction/blockvectorfunction/Makefile
     fem/operator/Makefile
     fem/operator/common/Makefile
     fem/operator/0order/Makefile
     fem/operator/1order/Makefile
     fem/operator/2order/Makefile
     fem/operator/matrix/Makefile
     fem/operator/projection/Makefile
     fem/io/Makefile
     fem/io/file/Makefile
     fem/io/visual/Makefile
     fem/io/visual/grape/Makefile
     fem/io/visual/grape/datadisp/Makefile
     fem/io/visual/grape/gridview/Makefile
     fem/io/visual/grape/functionview/Makefile
     fem/solver/Makefile 
     fem/solver/ode/Makefile 
     fem/solver/oemsolver/Makefile 
     dune-fem.pc])
AC_OUTPUT
