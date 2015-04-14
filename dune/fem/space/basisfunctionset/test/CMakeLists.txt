set(GRIDTYPE  YASPGRID)
set(GRIDDIM  2)
set(POLORDER  3)
set(DIMRANGE  5)
set(USE_VERTICAL_DOF_ALIGNMENT  1)

add_definitions( "-D${GRIDTYPE}" )
add_definitions( "-DGRIDDIM=${GRIDDIM}" ) 
add_definitions( "-DPOLORDER=${POLORDER}" )
add_definitions( "-DDIMRANGE=${DIMRANGE}" )
add_definitions( "-DUSE_VERTICAL_DOF_ALIGNMENT=${USE_VERTICAL_DOF_ALIGNMENT}" )

# copy data to build source to make tests work         
configure_file(1dgrid.dgf ${CMAKE_CURRENT_BINARY_DIR}/1dgrid.dgf COPYONLY)
configure_file(2dgrid.dgf ${CMAKE_CURRENT_BINARY_DIR}/2dgrid.dgf COPYONLY)
configure_file(3dgrid.dgf ${CMAKE_CURRENT_BINARY_DIR}/3dgrid.dgf COPYONLY)

add_executable(test_vectorialbasisfunctionset test-vectorialbasisfunctionset.cc)
dune_target_link_libraries(test_vectorialbasisfunctionset "${DUNE_LIBS};${LOCAL_LIBS}")
add_dune_mpi_flags(test_vectorialbasisfunctionset )
set_property(TARGET test_vectorialbasisfunctionset APPEND PROPERTY COMPILE_DEFINITIONS "USE_VERTICAL_DOF_ALIGNMENT=${USE_VERTICAL_DOF_ALIGNMENT}" )

add_executable(test_defaultbasisfunctionset test-defaultbasisfunctionset.cc)
dune_target_link_libraries(test_defaultbasisfunctionset "${DUNE_LIBS};${LOCAL_LIBS}")
add_dune_mpi_flags(test_defaultbasisfunctionset )
add_executable(test_simplebasisfunctionset test-simplebasisfunctionset.cc)
dune_target_link_libraries(test_simplebasisfunctionset "${DUNE_LIBS};${LOCAL_LIBS}")
add_dune_mpi_flags(test_simplebasisfunctionset )
add_executable(test_tuplebasisfunctionset test-tuplebasisfunctionset.cc)
dune_target_link_libraries(test_tuplebasisfunctionset "${DUNE_LIBS};${LOCAL_LIBS}")
add_dune_mpi_flags(test_tuplebasisfunctionset )

dune_add_test(test_vectorialbasisfunctionset test_defaultbasisfunctionset test_simplebasisfunctionset test_tuplebasisfunctionset)
dune_install(checkbasisfunctionset.hh)
