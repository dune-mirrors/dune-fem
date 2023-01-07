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
  testpetsc4py.hh
)

if( NOT DUNE_ENABLE_PYTHONBINDINGS )
  exclude_from_headercheck( ${HEADERS} )
endif()

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/fempy/py)
