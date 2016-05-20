add_subdirectory(grid)

set(HEADERS
  grid.hh
  space.hh
)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/fempy/py)
