if( NOT GRIDTYPE )
  set(GRIDTYPE YASPGRID)
endif()

if( NOT GRIDDIM )
  set(GRIDDIM 2)
endif()

dune_add_test( NAME quad_tester SOURCES main.cc cache_test.cc twist_test.cc point_test.cc
COMPILE_DEFINITIONS "${GRIDTYPE};GRIDDIM=${GRIDDIM}"
LINK_LIBRARIES dunefem )
