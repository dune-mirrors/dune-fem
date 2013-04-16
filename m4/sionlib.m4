# searches for sionlib-headers and libs

AC_DEFUN([DUNE_PATH_SIONLIB],[
  AC_REQUIRE([AC_PROG_CC])

  AC_ARG_WITH(sionlib,
    AC_HELP_STRING([--with-sionlib=PATH],[directory with SIONLIB inside]))
  AC_ARG_WITH(sionlibsuffix,
    AC_HELP_STRING([--with-sionlibsuffixe=SUFFIX],[suffix to the library name , e.g. gcc or something (default = "_64")]))
  AC_ARG_WITH(sionlib-includedir,
    AC_HELP_STRING([--with-sionlib-includedir=PATH],[directory with SIONLIB headers inside]))
  AC_ARG_WITH(sionlib-libdir,
    AC_HELP_STRING([--with-sionlib-libdir=PATH],[directory with SIONLIB libraries inside]))

# store old values
ac_save_LDFLAGS="$LDFLAGS"
ac_save_CPPFLAGS="$CPPFLAGS"
ac_save_LIBS="$LIBS"

SIONLIBYES=0
## do nothing if no --with-sionlib was supplied
if test x$with_sionlib != x && test x$with_sionlib != xno ; then
  SIONLIBYES=1
fi
SIONLIBSUFFIX="_64"
## do nothing if no --with-sionlibnsuffix was supplied
if test x$with_sionlibsuffix != x ; then
  SIONLIBSUFFIX=$with_sionlibsuffix
fi
if test x$with_sionlib_includedir != x && test x$with_sionlib_includedir != xno ; then
  SIONLIBYES=1
fi
if test x$with_sionlib_libdir != x && test x$with_sionlib_libdir != xno ; then
  SIONLIBYES=1
fi

if test x$SIONLIBYES = x1 ; then

  # is --with-sionlib=bla used?
  if test "x$with_sionlib" != x ; then
    SIONLIBROOT=`cd $with_sionlib && pwd`
	  if ! test -d $SIONLIBROOT;  then
      AC_MSG_WARN([SIONLIB directory $with_sionlib does not exist])
  	fi

    if test "x$SIONLIBROOT" = x; then
      # use some default value...
      SIONLIBROOT="/usr/local/sionlib"
    fi

    SIONLIB_LIB_PATH="$SIONLIBROOT/lib"
    SIONLIB_INCLUDE_PATH="$SIONLIBROOT/include"
  else 
    if test "x$with_sionlib_includedir" != x ; then 
      SIONLIB_INCLUDE_PATH=`cd $with_sionlib_includedir && pwd`
      if ! test -d $SIONLIB_INCLUDE_PATH;  then
        AC_MSG_WARN([SIONLIB directory $with_sionlib_includedir does not exist])
      fi
    fi
    if test "x$with_sionlib_libdir" != x ; then 
      SIONLIB_LIB_PATH=`cd $with_sionlib_libdir && pwd`
      if ! test -d $SIONLIB_LIB_PATH;  then
        AC_MSG_WARN([SIONLIB directory $with_sionlib_libdir does not exist])
      fi
    fi
  fi

  # set variables so that tests can use them
  REM_CPPFLAGS=$CPPFLAGS

  LDFLAGS="$LDFLAGS -L$SIONLIB_LIB_PATH"
  SIONLIB_INC_FLAG="-I$SIONLIB_INCLUDE_PATH -DENABLE_SIONLIB=1"
  CPPFLAGS="$CPPFLAGS $SIONLIB_INC_FLAG $MPI_CPPFLAGS"

  # check for header
  AC_LANG_PUSH([C])
  AC_CHECK_HEADERS([sion.h], 
     [SIONLIB_CPPFLAGS="$SIONLIB_INC_FLAG"
    HAVE_SIONLIB="1"],
    AC_MSG_WARN([sion.h not found in $SIONLIB_INCLUDE_PATH]))
   
  CPPFLAGS="$REM_CPPFLAGS"
  REM_CPPFLAGS=

  REM_LDFLAGS=$LDFLAGS

  # if header is found...
  if test x$HAVE_SIONLIB = x1 ; then
    LIBNAME="sion""$SIONLIBSUFFIX"
    AC_CHECK_LIB($LIBNAME,[main],
    [SIONLIB_LIBS="-l$LIBNAME"
           SIONLIB_LDFLAGS="-L$SIONLIB_LIB_PATH"],
	  [HAVE_SIONLIB="0"
	  AC_MSG_WARN(lib$LIBNAME not found!)])
  fi

  if test x$HAVE_SIONLIB = x1 ; then
    LIBNAME="sionser""$SIONLIBSUFFIX"
    AC_CHECK_LIB($LIBNAME,[main],
    [SIONLIB_LIBS="$SIONLIB_LIBS -l$LIBNAME"
     SIONLIB_LDFLAGS="-L$SIONLIB_LIB_PATH"],
	  [HAVE_SIONLIB="0"
	  AC_MSG_WARN(lib$LIBNAME not found!)])
  fi

  LDFLAGS=$REM_LDFLAGS
  AC_LANG_POP

## end of sionlib check (--without wasn't set)
fi

# survived all tests?
if test x$HAVE_SIONLIB = x1 ; then
  AC_SUBST(SIONLIB_LIBS, $SIONLIB_LIBS)
  AC_SUBST(SIONLIB_LDFLAGS, $SIONLIB_LDFLAGS)
  AC_SUBST(SIONLIB_CPPFLAGS, $SIONLIB_CPPFLAGS)
  AC_DEFINE(HAVE_SIONLIB, ENABLE_SIONLIB,
    [This is only true if sionlib-library was found by configure 
     _and_ if the application uses the SIONLIB_CPPFLAGS])

  # add to global list
  DUNE_ADD_ALL_PKG([SIONlib], [\${SIONLIB_CPPFLAGS}],
                       [\${SIONLIB_LDFLAGS}], [\${SIONLIB_LIBS}])

  # set variable for summary
  with_sionlib="yes"

else
  AC_SUBST(SIONLIB_LIBS, "")
  AC_SUBST(SIONLIB_LDFLAGS, "")
  AC_SUBST(SIONLIB_CPPFLAGS, "")

  # set variable for summary
  with_sionlib="no"
fi
  
# also tell automake
AM_CONDITIONAL(SIONLIB, test x$HAVE_SIONLIB = x1)

# reset old values
LIBS="$ac_save_LIBS"
CPPFLAGS="$ac_save_CPPFLAGS"
LDFLAGS="$ac_save_LDFLAGS"

  DUNE_ADD_SUMMARY_ENTRY([SIONLIB],[$with_sionlib])
])
