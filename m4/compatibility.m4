AC_DEFUN([DUNE_FEM_CHECK_COMPATIBILITY],[
  AC_ARG_ENABLE([compatibility],
    AC_HELP_STRING([--enable-compatibility],
                   [try to be compatible with older dune-fem versions]))
  AS_IF([ test "x$enable_compatibility" == "xyes" ],[
    AC_DEFINE([DUNE_FEM_COMPATIBILITY], [1], [Shall dune-fem be compatible with older versions?])
  ],[
    AC_DEFINE([DUNE_FEM_COMPATIBILITY], [0], [Shall dune-fem be compatible with older versions?])
  ])
])
