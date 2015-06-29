set( TESTS
  odetest
  test-timeprovider
)

foreach( test ${TESTS} )
  add_executable( ${test} EXCLUDE_FROM_ALL ${test}.cc )
  dune_target_link_libraries( ${test} "${DUNE_LIBS};${LOCAL_LIBS}" )
  add_test( ${test} ${test} )
endforeach()
