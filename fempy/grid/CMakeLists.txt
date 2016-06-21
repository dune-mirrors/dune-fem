set(HEADERS
  adaptation.hh
  adaptivedofvector.hh
  virtualizedrestrictprolong.hh
)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/fempy/grid)
