# some default flags needed very often
set( DEFAULTFLAGS "COUNT_FLOPS;YASPGRID;GRIDDIM=2" )

# dgl2projection tests for different discrete functions
set( tests eigen petsc )
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
