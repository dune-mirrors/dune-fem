set(HEADERS
  adaptation.hh
  discretefunctionmanager.hh
  virtualizedrestrictprolong.hh
)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/fempy/grid)
