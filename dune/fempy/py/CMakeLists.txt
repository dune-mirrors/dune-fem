add_subdirectory(common)
add_subdirectory(function)
add_subdirectory(grid)

set(HEADERS
  discretefunction.hh
  gridview.hh
  integrandsbase.hh
  integrands.hh
  operator.hh
  scheme.hh
  space.hh
  spaceadapt.hh
  ufllocalfunction.hh
)

if( NOT dune-python_FOUND )
  exclude_from_headercheck( ${HEADERS} )
endif()

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/fempy/py)
