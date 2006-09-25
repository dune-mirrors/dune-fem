AC_PREREQ(2.50)
AC_INIT(visual, 0.1, dune@mathematik.uni-freiburg.de)
AC_CONFIG_AUX_DIR([./])
AM_INIT_AUTOMAKE(visual, 0.1, dune@mathematik.uni-freiburg.de)


# check if we're in the right directory
AC_CONFIG_SRCDIR([grape/datadisp/datadisp.cc])

# which header-file is created
AM_CONFIG_HEADER(config.h)

# default install path is current directory 
AC_PREFIX_DEFAULT([`pwd`])

# we need no more than the standard DUNE-stuff
DUNE_CHECK_ALL_M([dune-common])
DUNE_CHECK_ALL_M([dune-grid])

if test x$HAVE_DUNE_COMMON != x1 ; then
  AC_MSG_ERROR([Can't work without the dune-common module. Maybe you have to supply your dune-common directory as --with-dunecommon=dir])
fi

if test x$HAVE_DUNE_GRID != x1 ; then
  AC_MSG_ERROR([Can't work without the dune-grid module. Maybe you have to supply your dune-common directory as --with-dunegrid=dir])
fi

# defines libtool to build the lib
AC_PROG_LIBTOOL

# implicitly set the Dune-flags everywhere
CPPFLAGS="$DUNE_CPPFLAGS"
LDFLAGS="$DUNE_LDFLAGS"
LIBS="$DUNE_LIBS"

# generate Makefiles
AC_CONFIG_FILES([Makefile
                 grape/Makefile
                 grape/gridview/Makefile
                 grape/datadisp/Makefile
                 grape/functionview/Makefile
               ])
AC_OUTPUT

# show search results
DUNE_SUMMARY_ALL
