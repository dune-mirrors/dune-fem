#! /bin/sh
# checks for dune-headers and everything they need

AC_DEFUN([DUNE_CHECK_OPTIONAL],[
  AC_REQUIRE([AC_PROG_CXX])
  AC_REQUIRE([AC_PROG_CXXCPP])
  AC_REQUIRE([AC_PROG_LIBTOOL])
  AC_REQUIRE([DUNE_DISABLE_LIBCHECK])

  # ____DUNE_CHECK_MODULES_____ ($1)

  m4_pushdef([_dune_name], [$1])
  m4_pushdef([_dune_module], [m4_translit(_dune_name, [-], [_])])
  m4_pushdef([_dune_header], [$2])
  m4_pushdef([_dune_ldpath], [lib])
  m4_pushdef([_dune_lib],    [m4_translit(_dune_name, [-], [])])
  m4_pushdef([_dune_symbol], [$3])
  m4_pushdef([_DUNE_MODULE], [m4_toupper(_dune_module)])

  # switch tests to c++
  AC_LANG_PUSH([C++])

  # the usual option...
  AC_ARG_WITH(_dune_name,
    AC_HELP_STRING([--with-_dune_name=PATH],[_dune_module directory]))

  # backup of flags
  ac_save_CPPFLAGS="$CPPFLAGS"
  ac_save_LIBS="$LIBS"
  ac_save_LDFLAGS="$LDFLAGS"
  CPPFLAGS=""
  LIBS=""

  ##
  ## Where is the module $1?
  ##

  AC_MSG_CHECKING([for $1 installation or source tree])

  # is a directory set?
  if test x$with_[]_dune_module = x ; then
    #
    # search module $1 via pkg-config
    #
    with_[]_dune_module="global installation"
    if test -z "$PKG_CONFIG"; then
      AC_MSG_RESULT([failed])
      AC_MSG_NOTICE([could not search for module _dune_name])
      AC_MSG_WARN([pkg-config is required for using installed modules])
    fi
    if AC_RUN_LOG([$PKG_CONFIG --exists --print-errors "$1"]); then
      _DUNE_MODULE[]_CPPFLAGS="`$PKG_CONFIG --cflags _dune_name`" 2>/dev/null
      _DUNE_MODULE[]_ROOT="`$PKG_CONFIG --variable=prefix _dune_name`" 2>/dev/null 
      ifelse(_dune_symbol,,,[
        _DUNE_MODULE[]_LDFLAGS="`$PKG_CONFIG --variable=libdir _dune_name`" 2>/dev/null 
        _DUNE_MODULE[]_LIBS="-l[]_dune_lib"
      ])
      dune_is_installed=1
      AC_MSG_RESULT([
        global installation in $_DUNE_MODULE[]_ROOT])
    else
      AC_MSG_RESULT([not found])
    fi
  else
    #
    # path for module $1 is specified via command line
    #
    if test -d $with_[]_dune_module ; then
      # expand tilde / other stuff
      _DUNE_MODULE[]_ROOT=`cd $with_[]_dune_module && pwd`

      # expand search path (otherwise empty CPPFLAGS)
      if test -d $_DUNE_MODULE[]_ROOT/include/dune; then
        # Dune was installed into directory given by with-dunecommon
        dune_is_installed=1
        _DUNE_MODULE[]_CPPFLAGS="-I$_DUNE_MODULE[]_ROOT/include"
      else
        _DUNE_MODULE[]_CPPFLAGS="-I$_DUNE_MODULE[]_ROOT"
      fi
      ifelse(_dune_symbol,,,[
        _DUNE_MODULE[]_LDFLAGS="-L$_DUNE_MODULE[]_ROOT/lib"
        _DUNE_MODULE[]_LIBS="-l[]_dune_lib"
      ])
      # set expanded module path
      with_[]_dune_module="$_DUNE_MODULE[]_ROOT"
      AC_MSG_RESULT([
        found in $_DUNE_MODULE[]_ROOT])
    else
      AC_MSG_RESULT([not found])
      AC_MSG_ERROR([_dune_name-directory $with_[]_dune_module does not exist])
    fi
  fi

  DUNE_CPPFLAGS="$DUNE_CPPFLAGS $_DUNE_MODULE[]_CPPFLAGS"
  CPPFLAGS="$DUNE_CPPFLAGS"
  SET_CPPFLAGS="$_DUNE_MODULE[]_CPPFLAGS"

  ##  
  ## check for an arbitrary header
  ##
  AC_CHECK_HEADER([dune/[]_dune_header],
    [HAVE_[]_DUNE_MODULE=1
     _DUNE_MODULE[]_CPPFLAGS="$SET_CPPFLAGS"],
    [HAVE_[]_DUNE_MODULE=0
     _DUNE_MODULE[]_CPPFLAGS=""
     AC_MSG_WARN([$with_[]_dune_module does not seem to contain a valid _dune_name (dune/[]_dune_header not found)])]
  )

  ##
  ## check for lib (if lib name was provided)
  ##
  ifelse(_dune_symbol,,AC_MSG_NOTICE([_dune_name does not provide libs]),[
    if test "x$enable_dunelibcheck" != "xyes"; then
      AC_MSG_WARN([library check for _dune_name is disabled. DANGEROUS!])
    fi
    if test x$HAVE_[]_DUNE_MODULE != x -a x$enable_dunelibcheck == "xyes"; then

      # save current LDFLAGS
      ac_save_LDFLAGS="$LDFLAGS"
      ac_save_LIBS="$LIBS"
      ac_save_CXX="$CXX"
      HAVE_[]_DUNE_MODULE=0

      # define LTCXXCOMPILE like it will be defined in the Makefile
      LTCXXLINK="$srcdir/libtool --tag=CXX --mode=link $CXX $CXXFLAGS $LDFLAGS"
      CXX="$LTCXXLINK"

      # use module LDFLAGS
      LDFLAGS="$LDFLAGS $DUNE_LDFLAGS $_DUNE_MODULE[]_LDFLAGS"
      LIBS="$DUNE_LIBS $_DUNE_MODULE[]_LIBS"

      AC_MSG_CHECKING([for lib[]_dune_lib])
      AC_TRY_LINK(dnl
        [#]include<dune/[]_dune_header>,
        _dune_symbol,
          [AC_MSG_RESULT([yes])
           HAVE_[]_DUNE_MODULE=1
           _DUNE_MODULE[]_LIBS="$LIBS"],
          [AC_MSG_RESULT([no])
           HAVE_[]_DUNE_MODULE=0
           AC_MSG_ERROR([$with_[]_dune_module does not seem to contain a valid _dune_name (failed to link with lib[]_dune_lib[].la)])]
      )

      # reset variables
      LDFLAGS="$ac_save_LDFLAGS"
      LIBS="$ac_save_LIBS"
      CXX="$ac_save_CXX"
    fi
  ])

  # did we succeed?
  if test x$HAVE_[]_DUNE_MODULE = x1 ; then
    # set variables for our modules
    AC_SUBST(_DUNE_MODULE[]_CPPFLAGS, "$_DUNE_MODULE[]_CPPFLAGS")
    AC_SUBST(_DUNE_MODULE[]_LDFLAGS, "$_DUNE_MODULE[]_LDFLAGS")
    AC_SUBST(_DUNE_MODULE[]_LIBS, "$_DUNE_MODULE[]_LIBS")
    AC_SUBST(_DUNE_MODULE[]_ROOT, "$_DUNE_MODULE[]_ROOT")
    AC_DEFINE(HAVE_[]_DUNE_MODULE, 1, [Define to 1 if _dune_module was found])

    # set DUNE_* variables
    AC_SUBST(DUNE_CPPFLAGS, "$DUNE_CPPFLAGS")
    AC_SUBST(DUNE_LDFLAGS, "$DUNE_LDFLAGS $_DUNE_MODULE[]_LDFLAGS")
    AC_SUBST(DUNE_LIBS, "$DUNE_LIBS $_DUNE_MODULE[]_LIBS")
    
    # add to global list
    # only add my flags other flags are added by other packages 
    DUNE_PKG_CPPFLAGS="$DUNE_PKG_CPPFLAGS $_DUNE_MODULE[]_CPPFLAGS"
    # also set all pkg flags here
    ALL_PKG_CPPFLAGS="$ALL_PKG_CPPFLAGS $_DUNE_MODULE[]_CPPFLAGS" 
    
    DUNE_PKG_LIBS="$DUNE_PKG_LIBS $LIBS"
    # also set all_pkg_libs here 
    ALL_PKG_LIBS="$ALL_PKG_LIBS $LIBS"
    
    DUNE_PKG_LDFLAGS="$DUNE_PKG_LDFLAGS $_DUNE_MODULE[]_LDFLAGS"
    ALL_PKG_LDFLAGS="$ALL_PKG_LDFLAGS $_DUNE_MODULE[]_LDFLAGS"

    with_[]_dune_module="yes"
  else
    with_[]_dune_module="no"
    AC_MSG_WARN([could not find required module _dune_name])
  fi

  # reset previous flags
  CPPFLAGS="$ac_save_CPPFLAGS"
  LIBS="$ac_save_LIBS"

  # remove local variables
  m4_popdef([_dune_name])
  m4_popdef([_dune_module])
  m4_popdef([_dune_header])
  m4_popdef([_dune_ldpath])
  m4_popdef([_dune_lib])
  m4_popdef([_dune_symbol])
  m4_popdef([_DUNE_MODULE])

  # restore previous language settings (leave C++)
  AC_LANG_POP([C++])
])
