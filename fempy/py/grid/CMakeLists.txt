set(HEADERS
  function.hh
  gridpart.hh
  hierarchical.hh
  numpy.hh
  restrictprolong.hh
)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/fempy/py/grid)
