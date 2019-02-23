# copy data to build source to make tests work         
configure_file(testparameterfile ${CMAKE_CURRENT_BINARY_DIR}/testparameterfile COPYONLY)
configure_file(param_DarcyStokes ${CMAKE_CURRENT_BINARY_DIR}/param_DarcyStokes COPYONLY)
configure_file(param_NavierStokes ${CMAKE_CURRENT_BINARY_DIR}/param_NavierStokes COPYONLY)
configure_file(bc.sh ${CMAKE_CURRENT_BINARY_DIR}/bc.sh COPYONLY)
configure_file(date.sh ${CMAKE_CURRENT_BINARY_DIR}/date.sh COPYONLY)

dune_add_test( NAME parametertest SOURCES parametertest.cc LINK_LIBRARIES dunefem )
dune_add_test( NAME parameterdicttest SOURCES parameterdicttest.cc LINK_LIBRARIES dunefem )
