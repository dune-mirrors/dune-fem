#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
# $Id$

# TODO:
#
# - put wrapper around mpi-check for consistency with other packages?

AC_PREREQ(2.50)
AC_INIT(fvcommon, 0.0, dune@mathematik.uni-freiburg.de)
AM_INIT_AUTOMAKE(fvcommon, 0.0, dune@mathematik.uni-freiburg.de)

# check if we're in the right directory
AC_CONFIG_SRCDIR([space/dgspace.hh])
# which header-file is created
AM_CONFIG_HEADER(config.h)

# we need no more than the standard DUNE-stuff
DUNE_CHECK_ALL_M([dunecommon])
DUNE_CHECK_ALL_M([dunegrid])

# implicitly set the Dune-flags everywhere
CPPFLAGS="$DUNE_CPPFLAGS"
LDFLAGS="$DUNE_LDFLAGS"
LIBS="$DUNE_LIBS"

# generate Makefiles
AC_CONFIG_FILES([Makefile
		 basefunctions/Makefile
		 basefunctions/test/Makefile
		 caching/Makefile
		 caching/test/Makefile
		 macrogridparser/Makefile
		 macrogridparser/test/Makefile
                 examples/Makefile
                 examples/poisson/Makefile
              	 examples/transport/Makefile
                 examples/vectorial/Makefile
                 examples/l2projection/Makefile
		 examples/localdg/Makefile
		 examples/localdg/advectdiffusion/Makefile
		 examples/pass/Makefile
		 pass/Makefile
		 pass/test/Makefile
		 space/Makefile
		 space/test/Makefile
		 quadrature/Makefile
		 quadrature/test/Makefile
               ])
AC_OUTPUT

# show search results
DUNE_SUMMARY_ALL
