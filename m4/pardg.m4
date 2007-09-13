# searches for pardg-headers and libs

AC_DEFUN([DUNE_PATH_PARDG],[
  AC_REQUIRE([AC_PROG_CXX])

  AC_ARG_WITH(pardg,
    AC_HELP_STRING([--with-pardg=PATH],[directory with ParDG inside]))
# do not use pardg debug lib 

# store old values
ac_save_LDFLAGS="$LDFLAGS"
ac_save_CPPFLAGS="$CPPFLAGS"
ac_save_LIBS="$LIBS"

## do nothing if no --with-pardg was supplied
if test x$with_pardg != x && test x$with_pardg != xno ; then

  if test x$with_pardg == xyes ; then
    AC_MSG_ERROR([You have to provide a directory --with-pardg=PATH])
  fi

  # is --with-pardg=bla used?
  if test "x$with_pardg" != x ; then
	if ! test -d $with_pardg; then
        AC_MSG_WARN([ParDG directory $with_pardg does not exist])
	else
    # expand tilde / other stuff
		PARDGROOT=`cd $with_pardg && pwd`
	fi
  fi

  if test "x$PARDGROOT" = x; then
    # use some default value...
    PARDGROOT="/usr/local/pardg"
  fi

  PARDG_LIB_PATH="$PARDGROOT/lib"
  PARDG_INCLUDE_PATH="$PARDGROOT/include"

  # set variables so that tests can use them
  REM_CPPFLAGS=$CPPFLAGS

  LDFLAGS="$LDFLAGS -L$PARDG_LIB_PATH"
  PARDG_INC_FLAG="-I$PARDG_INCLUDE_PATH -DENABLE_PARDG"
  CPPFLAGS="$CPPFLAGS $PARDG_INC_FLAG $MPI_CPPFLAGS"

  # check for header
  AC_LANG_PUSH([C++])
  AC_CHECK_HEADERS([function.hpp], 
     [PARDG_CPPFLAGS="$PARDG_INC_FLAG"
    HAVE_PARDG="1"],
    AC_MSG_WARN([function.hpp not found in $PARDG_INCLUDE_PATH]))
   
  CPPFLAGS="$REM_CPPFLAGS"
  REM_CPPFLAGS=

  REM_LDFLAGS=$LDFLAGS

  # if header is found...
  if test x$HAVE_PARDG = x1 ; then
    AC_CHECK_LIB(pardg,[main],
    [PARDG_LIBS="-lpardg"
           PARDG_LDFLAGS="-L$PARDG_LIB_PATH"
           LIBS="$LIBS $PARDG_LIBS"],
	  [HAVE_PARDG="0"
	  AC_MSG_WARN(libpardg not found!)])
  fi

  LDFLAGS=$REM_LDFLAGS
  AC_LANG_POP

## end of pardg check (--without wasn't set)
fi

# survived all tests?
if test x$HAVE_PARDG = x1 ; then
  AC_SUBST(PARDG_LIBS, $PARDG_LIBS)
  AC_SUBST(PARDG_LDFLAGS, $PARDG_LDFLAGS)
  AC_SUBST(PARDG_CPPFLAGS, $PARDG_CPPFLAGS)
  AC_DEFINE(HAVE_PARDG, ENABLE_PARDG,
    [This is only true if pardg-library was found by configure 
     _and_ if the application uses the PARDG_CPPFLAGS])

  # add to global list
  DUNE_PKG_LDFLAGS="$DUNE_PKG_LDFLAGS $PARDG_LDFLAGS"
  DUNE_PKG_LIBS="$DUNE_PKG_LIBS $PARDG_LIBS"
  DUNE_PKG_CPPFLAGS="$DUNE_PKG_CPPFLAGS $PARDG_CPPFLAGS"
  ALL_PKG_LDFLAGS="$ALL_PKG_LDFLAGS $PARDG_LDFLAGS"
  ALL_PKG_LIBS="$ALL_PKG_LIBS $PARDG_LIBS"
  ALL_PKG_CPPFLAGS="$ALL_PKG_CPPFLAGS $PARDG_CPPFLAGS"

  # set variable for summary
  with_pardg="yes"
else
  AC_SUBST(PARDG_LIBS, "")
  AC_SUBST(PARDG_LDFLAGS, "")
  AC_SUBST(PARDG_CPPFLAGS, "")

  # set variable for summary
  with_pardg="no"
fi
  
# also tell automake
AM_CONDITIONAL(PARDG, test x$HAVE_PARDG = x1)

# reset old values
LIBS="$ac_save_LIBS"
CPPFLAGS="$ac_save_CPPFLAGS"
LDFLAGS="$ac_save_LDFLAGS"

])
