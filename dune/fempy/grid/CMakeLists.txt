set(HEADERS
  adaptation.hh
  discretefunctionmanager.hh
  globalmapper.hh
  gridpartadapter.hh
  virtualizedrestrictprolong.hh
)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/fempy/grid)
