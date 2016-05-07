add_subdirectory(pybind11)

set(HEADERS
  gridfunction.hh
  grid.hh
)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/fempy)
