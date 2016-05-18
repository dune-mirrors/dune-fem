add_subdirectory(function)
add_subdirectory(pybind11)

set(HEADERS
  grid.hh
  parameter.hh
  pygrid.hh
  pygridfunction.hh
  pyspace.hh
  pyvtk.hh
)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/fempy)
