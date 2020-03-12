# set default grid type for make all

if( NOT GRIDTYPE )
  set(GRIDTYPE YASPGRID)
endif()

if( NOT GRIDDIM )
  set(GRIDDIM 3)
endif()


# copy data to build source to make tests work
set( TESTFILES 1dgrid.dgf 2dgrid.dgf 3dgrid.dgf  3dcubesimp.dgf parameter )
foreach( file ${TESTFILES} )
  configure_file(${file} ${CMAKE_CURRENT_BINARY_DIR}/${file} COPYONLY)
endforeach()


dune_add_test( NAME cachingquadrature SOURCES cachingquadrature.cc
COMPILE_DEFINITIONS "${GRIDTYPE};GRIDDIM=${GRIDDIM}"
LINK_LIBRARIES dunefem )

dune_add_test( NAME elementquadrature SOURCES elementquadrature.cc
COMPILE_DEFINITIONS "${GRIDTYPE};GRIDDIM=${GRIDDIM}"
LINK_LIBRARIES dunefem )

dune_add_test( NAME agglomerationquadrature SOURCES agglomerationquadrature.cc
COMPILE_DEFINITIONS "${GRIDTYPE};GRIDDIM=2"
LINK_LIBRARIES dunefem )
