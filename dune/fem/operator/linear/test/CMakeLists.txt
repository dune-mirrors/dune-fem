dune_install(checklinearoperator.hh)

# some default flags needed very often
set( DEFAULTFLAGS "COUNT_FLOPS;YASPGRID;GRIDDIM=2" )

dune_add_test( NAME test_matrices SOURCES test-matrices.cc COMPILE_DEFINITIONS "${DEFAULTFLAGS}"
LINK_LIBRARIES dunefem )

if( ${FEM_TORTURE_TESTS} )
  dune_add_test(
    SOURCES test-hierarchicallinearoperator.cc
    COMPILE_DEFINITIONS "YASPGRID" "GRIDDIM=2"
    LINK_LIBRARIES dunefem
    CMAKE_GUARD dune-istl_FOUND
  )
endif()
