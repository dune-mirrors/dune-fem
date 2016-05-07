set(DATAFILES
  dune-alugrid.db
  dune-grid.db
  dune-spgrid.db
)

foreach(DATAFILE ${DATAFILES})
  configure_file(${DATAFILE} ${DATAFILE} COPYONLY)
endforeach()
