set(DATAFILES
  femscheme.db
  navierscheme.db
)

foreach(DATAFILE ${DATAFILES})
  configure_file(${DATAFILE} ${DATAFILE} COPYONLY)
endforeach()
