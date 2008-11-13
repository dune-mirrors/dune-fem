AC_DEFUN([DUNE_FEM_CHECKS],[
  AC_REQUIRE([ACX_BLAS])
])

AC_DEFUN([DUNE_FEM_CHECK_MODULE],[
  DUNE_CHECK_MODULES([dune-fem], [fem/space/basefunctions/storageinterface.hh], [])
])
