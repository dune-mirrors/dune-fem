set(DATAFILES
  femscheme.db
)

foreach(DATAFILE ${DATAFILES})
  configure_file(${DATAFILE} ${DATAFILE} COPYONLY)
endforeach()
