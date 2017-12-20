set( HEADERS
  anisotropic.hh
  blockmapper.hh
  datahandle.hh
  default.hh
  legendre.hh
  localdofstorage.hh
  orthogonal.hh
  space.hh
)

install( FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/space/hpdg )
