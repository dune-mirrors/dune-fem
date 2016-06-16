set(DATAFILES
  femdiscretefunction.db
)

foreach(DATAFILE ${DATAFILES})
  configure_file(${DATAFILE} ${DATAFILE} COPYONLY)
endforeach()
