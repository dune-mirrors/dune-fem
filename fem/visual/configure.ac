AC_PREREQ(2.50)
AC_INIT(visual, 0.1, dune@mathematik.uni-freiburg.de)
AM_INIT_AUTOMAKE(visual, 0.1, dune@mathematik.uni-freiburg.de)

# check if we're in the right directory
AC_CONFIG_SRCDIR([grape/datadisp/datadisp.cc])
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
                 grape/Makefile
                 grape/gridview/Makefile
                 grape/datadisp/Makefile
               ])
AC_OUTPUT

# show search results
DUNE_SUMMARY_ALL
