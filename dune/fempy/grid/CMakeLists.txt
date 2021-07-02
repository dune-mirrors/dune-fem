set(HEADERS
  adaptation.hh
  discretefunctionmanager.hh
  globalmapper.hh
  virtualizedrestrictprolong.hh
)
exclude_from_headercheck( globalmapper.hh )

if( NOT DUNE_ENABLE_PYTHONBINDINGS )
  exclude_from_headercheck( ${HEADERS} )
endif()

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/fempy/grid)
