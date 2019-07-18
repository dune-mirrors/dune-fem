set(HEADERS
  gridfunctionview.hh
  simplegridfunction.hh
  virtualizedgridfunction.hh
)

if( NOT dune-python_FOUND )
  exclude_from_headercheck( ${HEADERS} )
endif()

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/fempy/function)
