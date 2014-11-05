set(TESTS  parametertest)

add_executable(parametertest parametertest.cc)
dune_target_link_libraries(parametertest "${DUNE_LIBS};${LOCAL_LIBS}")

# copy data to build source to make tests work         
configure_file(testparameterfile ${CMAKE_CURRENT_BINARY_DIR}/testparameterfile COPYONLY)
configure_file(param_DarcyStokes ${CMAKE_CURRENT_BINARY_DIR}/param_DarcyStokes COPYONLY)
configure_file(param_NavierStokes ${CMAKE_CURRENT_BINARY_DIR}/param_NavierStokes COPYONLY)
configure_file(bc.sh ${CMAKE_CURRENT_BINARY_DIR}/bc.sh COPYONLY)
configure_file(date.sh ${CMAKE_CURRENT_BINARY_DIR}/date.sh COPYONLY)

dune_add_test(i ${TESTS})
