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
AC_CONFIG_SRCDIR([tester/gridcheck/gridcheck.cc])
# which header-file is created
AM_CONFIG_HEADER(config.h)

# we need no more than the standard DUNE-stuff
DUNE_CHECK_ALL

# we can't possibly work without the dune headers...
if test x$HAVE_DUNE != x1 ; then
  AC_MSG_ERROR([Can't work without the DUNE-library. Maybe you'll have to pass a different path for --with-dune. See config.log for details])
fi

# implicitly set the Dune-flags everywhere
CPPFLAGS="$DUNE_CPPFLAGS"
LDFLAGS="$DUNE_LDFLAGS"
LIBS="$DUNE_LIBS"

# generate Makefiles
AC_CONFIG_FILES([Makefile
                 examples/Makefile
                 examples/poisson/Makefile
              	 examples/transport/Makefile
                 examples/vectorial/Makefile
                 examples/l2projection/Makefile
		 examples/localdg/Makefile
		 examples/localdg/advectdiffusion/Makefile
		 examples/pass/Makefile
		 examples/matrices/Makefile
                 tester/Makefile
                 tester/gridcheck/Makefile
                 tester/adaptation/Makefile
                 tester/gridcheck/unittests/Makefile
		 pass/Makefile
		 pass/test/Makefile
		 space/Makefile
		 quadrature/Makefile
		 quadrature/test/Makefile
               ])
AC_OUTPUT

# show search results
DUNE_SUMMARY_ALL

#set AC_PROG_CXX as desired by 
#tester/gridcheck/Makefile.am: The usual way to define `CXX' is to add `AC_PROG_CXX'
