# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.50)
DUNE_AC_INIT
AM_INIT_AUTOMAKE([foreign 1.9 tar-pax dist-bzip2])

# check if we're in the right directory
AC_CONFIG_SRCDIR([dune/fem/space/finitevolume.hh])
# which header-file is created
AM_CONFIG_HEADER([config.h])

# check all dune-module stuff
DUNE_CHECK_ALL

DUNE_SUMMARY_ALL

# preset variable to path such that #include <dune/...> works
AC_SUBST(AM_CPPFLAGS, '$(DUNE_CPPFLAGS) -I$(top_srcdir)')
AC_SUBST(AM_LDFLAGS, '$(DUNE_LDFLAGS) $(DUNE_LIBS)')
AC_SUBST([LOCAL_LIBS], '$(top_builddir)/dune/fem/libfem.la')

# implicitly set the Dune-flags everywhere
CPPFLAGS="$DUNE_CPPFLAGS"
LDFLAGS="$DUNE_LDFLAGS"
LIBS="$DUNE_LIBS $BLAS_LIBS"

AC_CONFIG_FILES([Makefile
     doc/Makefile
     doc/doxygen/Makefile
     doc/doxygen/pages/Makefile
     dune/Makefile
     dune/fem/Makefile
     dune/fem/pass/Makefile
     dune/fem/pass/test/Makefile
     dune/fem/pass/test/dgellipt/Makefile
     dune/fem/space/Makefile
     dune/fem/space/test/Makefile
     dune/fem/space/basefunctions/Makefile
     dune/fem/space/basisfunctionset/Makefile
     dune/fem/space/shapefunctionset/Makefile
     dune/fem/space/common/Makefile
     dune/fem/space/dofmapper/Makefile
     dune/fem/space/discontinuousgalerkin/Makefile
     dune/fem/space/discontinuousgalerkin/test/Makefile
     dune/fem/space/finitevolume/Makefile
     dune/fem/space/finitevolume/test/Makefile
     dune/fem/space/lagrange/Makefile
     dune/fem/space/lagrange/test/Makefile
     dune/fem/space/rannacherturek/Makefile
     dune/fem/space/rannacherturek/test/Makefile
     dune/fem/space/reducedbasisspace/Makefile
     dune/fem/space/reducedbasisspace/test/Makefile
     dune/fem/space/combinedspace/Makefile
     dune/fem/space/padaptivespace/Makefile
     dune/fem/space/mapper/Makefile
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
     dune/fem/function/blockvectors/Makefile
     dune/fem/function/blockvectordiscretefunction/Makefile
     dune/fem/function/petscdiscretefunction/Makefile
     dune/fem/operator/Makefile
     dune/fem/operator/common/Makefile
     dune/fem/operator/common/test/Makefile
     dune/fem/operator/0order/Makefile
     dune/fem/operator/1order/Makefile
     dune/fem/operator/2order/Makefile
     dune/fem/operator/matrix/Makefile
     dune/fem/operator/model/Makefile
     dune/fem/operator/projection/Makefile
     dune/fem/petsc/tests/Makefile
     dune/fem/petsc/tests/l2projection/Makefile
     dune/fem/petsc/experiments/petsc_matrix/Makefile
     dune/fem/petsc/experiments/Makefile
     dune/fem/petsc/experiments/commtest/Makefile
     dune/fem/petsc/experiments/parallel_laplace/Makefile
     dune/fem/petsc/Makefile
     dune/fem/io/Makefile
     dune/fem/io/test/Makefile
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
     dune/fem/misc/petsc/Makefile
     dune/fem/storage/Makefile
     dune/fem/test/Makefile
     dune/fem/gridpart/Makefile
     dune/fem/gridpart/common/Makefile
     dune/fem/gridpart/filter/Makefile
     dune/fem/gridpart/filteredgridpart/Makefile
     dune/fem/gridpart/geogridpart/Makefile
     dune/fem/gridpart/idgridpart/Makefile
     dune/fem/gridpart/test/Makefile
     lib/Makefile
     m4/Makefile
     make.inc
     dune-fem.pc])
AC_OUTPUT
