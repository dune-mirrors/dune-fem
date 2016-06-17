add_subdirectory(function)
add_subdirectory(py)
add_subdirectory(pybind11)
add_subdirectory(quadrature)

set(HEADERS
  grid.hh
  parameter.hh
)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/fempy)
