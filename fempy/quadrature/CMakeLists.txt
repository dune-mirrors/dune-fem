set(HEADERS
  cachingpoint.hh
  elementpoint.hh
  fempyquadratures.hh
)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/fempy/quadrature)
