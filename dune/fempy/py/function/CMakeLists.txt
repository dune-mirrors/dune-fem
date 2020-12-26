set(HEADERS
  grid.hh
)

if( NOT DUNE_ENABLE_PYTHONBINDINGS )
  exclude_from_headercheck( ${HEADERS} )
endif()

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/fempy/py/function)
