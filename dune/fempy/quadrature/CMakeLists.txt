set(HEADERS
  cachingpoint.hh
  elementpoint.hh
  fempyquadratures.hh
)

if( NOT dune-python_FOUND )
  exclude_from_headercheck( ${HEADERS} )
endif()

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/fempy/quadrature)
