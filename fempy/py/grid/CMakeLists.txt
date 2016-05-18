set(HEADERS
  entity.hh
  geometry.hh
  hierarchical.hh
)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/fempy/py/grid)
