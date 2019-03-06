set(HEADERS
  default.hh
  doerfler.hh
  localerror.hh
  maximum.hh
)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/marking)
