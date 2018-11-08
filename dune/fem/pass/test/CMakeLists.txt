if( ${TORTURE_TESTS} )
  dune_add_test( NAME test_insertoperatorpass SOURCES test-insertoperatorpass.cc LINK_LIBRARIES dunefem )
endif()