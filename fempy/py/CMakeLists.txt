add_subdirectory(common)
add_subdirectory(function)
add_subdirectory(grid)

set(HEADERS
  grid.hh
  space.hh
  discretefunction.hh
)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/fempy/py)
