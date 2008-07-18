# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.50)
DUNE_AC_INIT
AM_INIT_AUTOMAKE([foreign 1.9 tar-pax dist-bzip2])

# check if we're in the right directory
AC_CONFIG_SRCDIR([fem/space/dgspace/dgspace.hh])
# which header-file is created
AM_CONFIG_HEADER([config.h])

# we need no more than the standard DUNE-stuff
DUNE_CHECK_ALL_M

# parse dune-fem version
DUNE_FEM_PARSE_VERSION

# check dennis package 
DUNE_PATH_PARDG

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
AC_SUBST([LOCAL_LIBS], '$(top_builddir)/fem/libfem.la')

# implicitly set the Dune-flags everywhere
CPPFLAGS="$DUNE_CPPFLAGS"
LDFLAGS="$DUNE_LDFLAGS"
LIBS="$DUNE_LIBS $BLAS_LIBS"

# generate Makefiles
AC_CHECK_FILE([fem/devel-suite/Makefile.am],
              [AC_CONFIG_FILES([fem/devel-suite/Makefile
                                fem/devel-suite/poisson/Makefile
                                fem/devel-suite/elliptic/Makefile
                                fem/devel-suite/transport/Makefile
                                fem/devel-suite/vectorial/Makefile
                                fem/devel-suite/l2projection/Makefile
                                fem/devel-suite/localdg/Makefile
                                fem/devel-suite/localdg/advectdiffusion/Makefile
                                fem/devel-suite/pass/Makefile
                                fem/devel-suite/laplace/Makefile])],
              [])

AC_CHECK_FILE([fem/pass/test/Makefile.am],
              [AC_CONFIG_FILES([fem/pass/test/Makefile
                                fem/pass/test/dgellipt/Makefile])],
              [])
              
AC_CONFIG_FILES([Makefile
     doc/Makefile
     doc/doxygen/Makefile
     fem/Makefile
     fem/pass/Makefile
     fem/space/Makefile
     fem/space/test/Makefile
     fem/space/basefunctions/Makefile
     fem/space/dgspace/Makefile
     fem/space/fvspace/Makefile
     fem/space/lagrangespace/Makefile
     fem/space/lagrangespace/test/Makefile
     fem/space/reducedbasisspace/Makefile
     fem/space/reducedbasisspace/test/Makefile
     fem/space/combinedspace/Makefile
     fem/space/common/Makefile
     fem/quadrature/Makefile
     fem/quadrature/test/Makefile
     fem/quadrature/caching/Makefile
     fem/quadrature/caching/test/Makefile
     fem/function/Makefile
     fem/function/common/Makefile
     fem/function/localfunction/Makefile
     fem/function/adaptivefunction/Makefile
     fem/function/attachedfunction/Makefile
     fem/function/blockvectorfunction/Makefile
     fem/function/combinedfunction/Makefile
     fem/function/productfunction/Makefile
     fem/function/productfunction/test/Makefile
     fem/function/vectorfunction/Makefile
     fem/operator/Makefile
     fem/operator/common/Makefile
     fem/operator/common/test/Makefile
     fem/operator/0order/Makefile
     fem/operator/1order/Makefile
     fem/operator/2order/Makefile
     fem/operator/matrix/Makefile
     fem/operator/model/Makefile
     fem/operator/projection/Makefile
     fem/io/Makefile
     fem/io/file/Makefile
     fem/io/file/test/Makefile
     fem/io/streams/Makefile
     fem/io/streams/test/Makefile
     fem/io/matlab/Makefile
     fem/io/matlab/m/Makefile
     fem/io/matlab/test/Makefile
     fem/io/visual/Makefile
     fem/io/visual/grape/Makefile
     fem/io/visual/grape/datadisp/Makefile
     fem/io/visual/grape/gridview/Makefile
     fem/io/visual/grape/functionview/Makefile
     fem/solver/Makefile 
     fem/solver/test/Makefile 
     fem/solver/ode/Makefile 
     fem/solver/oemsolver/Makefile 
     fem/misc/Makefile
     fem/storage/Makefile
     fem/test/Makefile
     fem/gridpart/Makefile
     lib/Makefile
     scripts/check-header.sh
     fem/version.hh
     dune-fem.pc])
AC_OUTPUT
