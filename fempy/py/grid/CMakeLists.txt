set(HEADERS
  adaptation.hh
  function.hh
  gridpart.hh
  numpy.hh
  restrictprolong.hh
)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/fempy/py/grid)
