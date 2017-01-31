set(HEADERS
  sharedgeometry.hh
  unimplementedgeometry.hh
)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/gridpart/common)
