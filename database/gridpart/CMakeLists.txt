set(DATAFILES
  dune-fem.db
)

foreach(DATAFILE ${DATAFILES})
  configure_file(${DATAFILE} ${DATAFILE} COPYONLY)
endforeach()
