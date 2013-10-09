# searches for papi-headers and libs

AC_DEFUN([DUNE_PATH_PAPI],[
  AC_REQUIRE([AC_PROG_CC])

  AC_ARG_WITH(papi,
    AC_HELP_STRING([--with-papi=PATH],[directory with PAPI inside]))

# store old values
ac_save_LDFLAGS="$LDFLAGS"
ac_save_CPPFLAGS="$CPPFLAGS"
ac_save_LIBS="$LIBS"

PAPIYES=0
## do nothing if no --with-papi was supplied
if test x$with_papi != x && test x$with_papi != xno ; then
  PAPIYES=1
fi

if test x$PAPIYES = x1 ; then

  # is --with-papi=bla used?
  if test "x$with_papi" != x ; then
    PAPIROOT=`cd $with_papi && pwd`
	  if ! test -d $PAPIROOT;  then
      AC_MSG_WARN([PAPI directory $with_papi does not exist])
  	fi

    if test "x$PAPIROOT" = x; then
      # use some default value...
      PAPIROOT="/usr/local/papi"
    fi

    PAPI_LIB_PATH="$PAPIROOT/lib"
    PAPI_INCLUDE_PATH="$PAPIROOT/include"
  fi

  # set variables so that tests can use them
  REM_CPPFLAGS=$CPPFLAGS

  LDFLAGS="$LDFLAGS -L$PAPI_LIB_PATH"
  PAPI_INC_FLAG="-I$PAPI_INCLUDE_PATH -DENABLE_PAPI=1"
  CPPFLAGS="$CPPFLAGS $PAPI_INC_FLAG"

  # check for header
  AC_LANG_PUSH([C])
  AC_CHECK_HEADERS([papi.h], 
     [PAPI_CPPFLAGS="$PAPI_INC_FLAG"
    HAVE_PAPI="1"],
    AC_MSG_WARN([papi.h not found in $PAPI_INCLUDE_PATH]))
   
  CPPFLAGS="$REM_CPPFLAGS"
  REM_CPPFLAGS=

  REM_LDFLAGS=$LDFLAGS

  # if header is found...
  if test x$HAVE_PAPI = x1 ; then
    LIBNAME="papi"
    AC_CHECK_LIB($LIBNAME,[PAPI_flops],
    [PAPI_LIBS="-l$LIBNAME"
           PAPI_LDFLAGS="-L$PAPI_LIB_PATH"],
	  [HAVE_PAPI="0"
	  AC_MSG_WARN(lib$LIBNAME not found!)])
  fi

  LDFLAGS=$REM_LDFLAGS
  AC_LANG_POP
fi

# survived all tests?
if test x$HAVE_PAPI = x1 ; then
  AC_SUBST(PAPI_LIBS, $PAPI_LIBS)
  AC_SUBST(PAPI_LDFLAGS, $PAPI_LDFLAGS)
  AC_SUBST(PAPI_CPPFLAGS, $PAPI_CPPFLAGS)
  AC_DEFINE(HAVE_PAPI, ENABLE_PAPI,
    [This is only true if papi-library was found by configure 
     _and_ if the application uses the PAPI_CPPFLAGS])

  # add to global list
  DUNE_ADD_ALL_PKG([PAPI], [\${PAPI_CPPFLAGS}],
                   [\${PAPI_LDFLAGS}], [\${PAPI_LIBS}])

  # set variable for summary
  with_papi="yes"

else
  AC_SUBST(PAPI_LIBS, "")
  AC_SUBST(PAPI_LDFLAGS, "")
  AC_SUBST(PAPI_CPPFLAGS, "")

  # set variable for summary
  with_papi="no"
fi
  
# also tell automake
AM_CONDITIONAL(PAPI, test x$HAVE_PAPI = x1)

# reset old values
LIBS="$ac_save_LIBS"
CPPFLAGS="$ac_save_CPPFLAGS"
LDFLAGS="$ac_save_LDFLAGS"

DUNE_ADD_SUMMARY_ENTRY([PAPI],[$with_papi])
])
