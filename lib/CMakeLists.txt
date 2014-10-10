dune_add_library(dunefem 
    _DUNE_TARGET_OBJECTS:io_
    _DUNE_TARGET_OBJECTS:misc_ 
    _DUNE_TARGET_OBJECTS:space_ 
    _DUNE_TARGET_OBJECTS:solver_ 
    _DUNE_TARGET_OBJECTS:quadrature_
    ADD_LIBS dunegeometry)
add_dune_mpi_flags(OBJECT dunefem)
