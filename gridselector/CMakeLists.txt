set( GRIDDEFINITIONS
  dune-alugrid
  dune-grid
  dune-spgrid
)

install( FILES ${GRIDDEFINITIONS} DESTINATION ${CMAKE_INSTALL_DATADIR}/gridselector )
