dune_add_library(dunefem
    _DUNE_TARGET_OBJECTS:io_
    _DUNE_TARGET_OBJECTS:misc_
    _DUNE_TARGET_OBJECTS:space_
    _DUNE_TARGET_OBJECTS:solver_
    _DUNE_TARGET_OBJECTS:storage_
    _DUNE_TARGET_OBJECTS:quadrature_
    ADD_LIBS "${PAPI_LIBRARIES};${SIONLIB_LIBRARIES};${DUNE_LIBS}")

  ##does the following work?
  add_dune_mpi_flags(OBJECT dunefem)

install(FILES DESTINATION ${CMAKE_INSTALL_BINDIR}/../lib)
