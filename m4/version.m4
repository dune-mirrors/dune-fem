AC_DEFUN([DUNE_FEM_PARSE_VERSION],[
  DUNE_FEM_VERSION="$PACKAGE_VERSION"
  DUNE_FEM_VERSION_MAJOR="`echo $DUNE_FEM_VERSION | awk '{ split($[]0,s,"."); match(s[[1]],/[[[[:digit:]]]]*/,n); print n[[0]] }'`"
  DUNE_FEM_VERSION_MINOR="`echo $DUNE_FEM_VERSION | awk '{ split($[]0,s,"."); match(s[[2]],/[[[[:digit:]]]]*/,n); print n[[0]] }'`"
  DUNE_FEM_VERSION_REVISION="`echo $DUNE_FEM_VERSION | awk '{ split($[]0,s,"."); match(s[[3]],/[[[[:digit:]]]]*/,n); print n[[0]] }'`"

  if test x$DUNE_FEM_VERSION_MAJOR = x -o x$DUNE_FEM_VERSION_MINOR = x -o DUNE_FEM_VERSION_REVISION = x ; then
    AC_MSG_ERROR([Cannot parse dune-fem version: $DUNE_FEM_VERSION])
  fi
  AC_MSG_NOTICE([Parsed dune-fem version: $DUNE_FEM_VERSION_MAJOR.$DUNE_FEM_VERSION_MINOR.$DUNE_FEM_VERSION_REVISION])
  AC_DEFINE_UNQUOTED(DUNE_FEM_VERSION_MAJOR, $DUNE_FEM_VERSION_MAJOR, [Define to the major version of dune-fem])
  AC_DEFINE_UNQUOTED(DUNE_FEM_VERSION_MINOR, $DUNE_FEM_VERSION_MINOR, [Define to the minor version of dune-fem])
  AC_DEFINE_UNQUOTED(DUNE_FEM_VERSION_REVISION, $DUNE_FEM_VERSION_REVISION, [Define to the revision of dune-fem])

  DUNE_GRID_VERSION_MAJOR="`echo $DUNE_GRID_VERSION | awk '{ split($[]0,s,"."); match(s[[1]],/[[[[:digit:]]]]*/,n); print n[[0]] }'`"
  DUNE_GRID_VERSION_MINOR="`echo $DUNE_GRID_VERSION | awk '{ split($[]0,s,"."); match(s[[2]],/[[[[:digit:]]]]*/,n); print n[[0]] }'`"
  if test x$DUNE_GRID_VERSION_MAJOR = x -o DUNE_GRID_VERSION_MINOR = x ; then
    AC_MSG_ERROR([Cannot parse dune-grid version: $DUNE_GRID_VERSION])
  fi
  AC_MSG_NOTICE([Parsed dune-grid version: $DUNE_GRID_VERSION_MAJOR.$DUNE_GRID_VERSION_MINOR])
  if test \( $DUNE_GRID_VERSION_MAJOR -lt 1 \) -o \( \( $DUNE_GRID_VERSION_MAJOR -eq 1 \) -a \( $DUNE_GRID_VERSION_MINOR -lt 1 \) \) ; then
    AC_MSG_ERROR([Your 'dune-grid' version too old; at least version 1.1 is required.])
  fi
  AC_DEFINE_UNQUOTED(DUNE_GRID_VERSION_MAJOR, $DUNE_GRID_VERSION_MAJOR, [Define to the major version of dune-grid])
  AC_DEFINE_UNQUOTED(DUNE_GRID_VERSION_MINOR, $DUNE_GRID_VERSION_MINOR, [Define to the minor version of dune-grid])
  AC_DEFINE_UNQUOTED(DUNE_GRID_VERSION_REVISION, 0, [Define to the revision of dune-grid])
])
