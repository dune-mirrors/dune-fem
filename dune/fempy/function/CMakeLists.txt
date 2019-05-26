set(HEADERS
  gridfunctionview.hh
  simplegridfunction.hh
  virtualizedgridfunction.hh
)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/fempy/function)
