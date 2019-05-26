set(HEADERS
  compositegeometry.hh
  localfunctiongeometry.hh
  sharedgeometry.hh
  simplegeometry.hh
  unimplementedgeometry.hh
)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/gridpart/common)
