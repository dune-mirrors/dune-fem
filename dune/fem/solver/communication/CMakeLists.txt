set(HEADERS
  fem.hh
  hierarchical.hh
  owneroverlapcopy.hh
)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/solver/communication)
