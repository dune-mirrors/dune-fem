set( TESTS test-hybrid test-staticints)

#finally add all pkg flags and local libs to tests
foreach( test ${TESTS} )
  dune_add_test( NAME ${test} SOURCES ${test}.cc )
  target_link_dune_default_libraries( ${test} )
  dune_target_link_libraries( ${test} dunefem )
endforeach()
