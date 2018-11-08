dune_install(checklinearoperator.hh)

# some default flags needed very often
set( DEFAULTFLAGS "COUNT_FLOPS;YASPGRID;GRIDDIM=2" )

# dgl2projection tests for different discrete functions
set( tests petsc )
if( ${TORTURE_TESTS} )
  set( tests "${tests} eigen" )
endif()
foreach( matrix ${tests} )
  string( TOUPPER ${matrix} MATRIX )
  if( HAVE_${MATRIX} )
    dune_add_test( NAME test_${matrix}matrix SOURCES test-matrices.cc COMPILE_DEFINITIONS "${DEFAULTFLAGS};USE_${MATRIX}"
    LINK_LIBRARIES dunefem )
  endif()
endforeach()

dune_add_test( NAME test_spmatrix SOURCES test-matrices.cc COMPILE_DEFINITIONS "${DEFAULTFLAGS}"
LINK_LIBRARIES dunefem )

dune_add_test( NAME test_istlmatrix SOURCES test-matrices.cc COMPILE_DEFINITIONS "${DEFAULTFLAGS};USE_ISTL"
LINK_LIBRARIES dunefem )

if( ${TORTURE_TESTS} )
  dune_add_test(
    SOURCES test-hierarchicallinearoperator.cc
    COMPILE_DEFINITIONS "YASPGRID" "GRIDDIM=2"
    LINK_LIBRARIES dunefem
    CMAKE_GUARD dune-istl_FOUND
  )
endif()
