set(HEADERS
  entity.hh
  function.hh
  gridpart.hh
  geometry.hh
  hierarchical.hh
  numpy.hh
  restrictprolong.hh
  vtk.hh
)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/fempy/py/grid)
