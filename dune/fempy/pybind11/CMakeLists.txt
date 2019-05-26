set(HEADERS
  gridfunction.hh
  space.hh
  pybind11.hh
)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/fempy/pybind11)
