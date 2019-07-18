set(HEADERS
  adaptation.hh
  discretefunctionmanager.hh
  globalmapper.hh
  gridpartadapter.hh
  virtualizedrestrictprolong.hh
)

if( NOT dune-python_FOUND )
  exclude_from_headercheck( ${HEADERS} )
endif()

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/fempy/grid)
