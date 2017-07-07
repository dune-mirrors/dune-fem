# copy data to build source to make tests work
configure_file(parameter ${CMAKE_CURRENT_BINARY_DIR}/parameter COPYONLY)

dune_add_test( NAME odetest SOURCES odetest.cc LINK_LIBRARIES dunefem )

dune_add_test( NAME test_timeprovider SOURCES test-timeprovider.cc LINK_LIBRARIES dunefem )

dune_add_test( NAME newtontest SOURCES newtontest.cc LINK_LIBRARIES dunefem )

dune_add_test( NAME inverseoperatortest SOURCES inverseoperatortest.cc LINK_LIBRARIES dunefem )
