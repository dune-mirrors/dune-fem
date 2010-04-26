# Define DUNE_OUTPUT_GRIDTYPE, if dune-grid does not already define it
AC_DEFUN([DUNE_OUTPUT_GRIDTYPE],[
  AC_DEFUN([DUNE_OUTOUT_GRIDTYPE],[
    GRIDTYPE_HH="gridtype.hh"
    AC_MSG_NOTICE([creating compatibility version of $GRIDTYPE_HH])

    echo "#include <dune/grid/utility/gridtype.hh> > $GRIDTYPE_HH

    DGFGRIDTYPE_HH="dgfgridtype.hh"
    AC_MSG_NOTICE([creating compatibility version of $DGFGRIDTYPE_HH])

    echo "#include <dune/grid/io/file/dgfparser/dgfparser/dgfgridtype.hh> > $DGFGRIDTYPE_HH
  ])
])

AC_DEFUN([DUNE_FEM_CHECKS],[
  AC_REQUIRE([ACX_BLAS])
  AC_REQUIRE([DUNE_PATH_PARDG])
  AC_REQUIRE([DUNE_PATH_UMFPACK])
])

AC_DEFUN([DUNE_FEM_CHECK_MODULE],[
  DUNE_CHECK_MODULES([dune-fem], [fem/space/dgspace/orthonormalbase_mod.hh], [OrthonormalBase_1D::eval_line])
])
