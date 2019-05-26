add_subdirectory(common)
add_subdirectory(function)
add_subdirectory(grid)

set(HEADERS
  discretefunction.hh
  gridview.hh
  integrands.hh
  operator.hh
  scheme.hh
  space.hh
  spaceadapt.hh
  ufllocalfunction.hh
)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/fempy/py)
