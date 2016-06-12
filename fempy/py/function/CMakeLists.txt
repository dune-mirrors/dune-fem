set(HEADERS
  discrete.hh
  grid.hh
)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/fempy/py/function)
