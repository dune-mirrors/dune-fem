# Process this file with autoconf to produce a configure script.

# Stupid comment line

AC_PREREQ(2.50)
DUNE_AC_INIT
AM_INIT_AUTOMAKE([foreign 1.9 tar-pax dist-bzip2])

# check if we're in the right directory
AC_CONFIG_SRCDIR([dune/fem/space/dgspace/dgspace.hh])
# which header-file is created
AM_CONFIG_HEADER([config.h])

# check all dune-module stuff
DUNE_CHECK_ALL

AC_ARG_ENABLE(
  compatibility,
  AC_HELP_STRING([--enable-compatibility],
                 [try to be compatible with older dune-fem versions]))
if test x$enable_compatibility == xyes ; then
  AC_DEFINE(DUNE_FEM_COMPATIBILITY,1,[Shall dune-fem be compatible with older versions?])
else
  AC_DEFINE(DUNE_FEM_COMPATIBILITY,0,[Shall dune-fem be compatible with older versions?])
fi


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
AC_SUBST([LOCAL_LIBS], '$(top_builddir)/dune/fem/libfem.la')

# implicitly set the Dune-flags everywhere
CPPFLAGS="$DUNE_CPPFLAGS"
LDFLAGS="$DUNE_LDFLAGS"
LIBS="$DUNE_LIBS $BLAS_LIBS"

AC_CHECK_FILE([dune/fem/pass/test/Makefile.am],
              [AC_CONFIG_FILES([dune/fem/pass/test/Makefile
                                dune/fem/pass/test/dgellipt/Makefile])],
              [])
              
AC_CONFIG_FILES([Makefile
     doc/Makefile
     doc/doxygen/Makefile
     doc/doxygen/pages/Makefile
     dune/Makefile
     dune/fem/Makefile
     dune/fem/pass/Makefile
     dune/fem/space/Makefile
     dune/fem/space/test/Makefile
     dune/fem/space/basefunctions/Makefile
     dune/fem/space/dgspace/Makefile
     dune/fem/space/fvspace/Makefile
     dune/fem/space/lagrangespace/Makefile
     dune/fem/space/lagrangespace/test/Makefile
     dune/fem/space/pqlagrangespace/Makefile
     dune/fem/space/pqlagrangespace/test/Makefile
     dune/fem/space/reducedbasisspace/Makefile
     dune/fem/space/reducedbasisspace/test/Makefile
     dune/fem/space/combinedspace/Makefile
     dune/fem/space/common/Makefile
     dune/fem/quadrature/Makefile
     dune/fem/quadrature/test/Makefile
     dune/fem/quadrature/caching/Makefile
     dune/fem/quadrature/caching/test/Makefile
     dune/fem/function/Makefile
     dune/fem/function/common/Makefile
     dune/fem/function/localfunction/Makefile
     dune/fem/function/adaptivefunction/Makefile
     dune/fem/function/attachedfunction/Makefile
     dune/fem/function/blockvectorfunction/Makefile
     dune/fem/function/combinedfunction/Makefile
     dune/fem/function/productfunction/Makefile
     dune/fem/function/productfunction/test/Makefile
     dune/fem/function/vectorfunction/Makefile
     dune/fem/operator/Makefile
     dune/fem/operator/common/Makefile
     dune/fem/operator/common/test/Makefile
     dune/fem/operator/0order/Makefile
     dune/fem/operator/1order/Makefile
     dune/fem/operator/2order/Makefile
     dune/fem/operator/matrix/Makefile
     dune/fem/operator/model/Makefile
     dune/fem/operator/projection/Makefile
     dune/fem/io/Makefile
     dune/fem/io/file/Makefile
     dune/fem/io/file/test/Makefile
     dune/fem/io/streams/Makefile
     dune/fem/io/streams/test/Makefile
     dune/fem/io/matlab/Makefile
     dune/fem/io/matlab/m/Makefile
     dune/fem/io/matlab/test/Makefile
     dune/fem/io/visual/Makefile
     dune/fem/io/visual/grape/Makefile
     dune/fem/io/visual/grape/datadisp/Makefile
     dune/fem/solver/Makefile 
     dune/fem/solver/test/Makefile 
     dune/fem/solver/ode/Makefile 
     dune/fem/solver/oemsolver/Makefile 
     dune/fem/misc/Makefile
     dune/fem/storage/Makefile
     dune/fem/test/Makefile
     dune/fem/gridpart/Makefile
     lib/Makefile
     m4/Makefile
     make.inc
     dune-fem.pc])
AC_OUTPUT
DUNE_OUTPUT_GRIDTYPE
