add_subdirectory(grid)

set(HEADERS
  grid.hh
  gridfunction.hh
  space.hh
  vtk.hh
)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/fempy/py)
