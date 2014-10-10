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


set(TESTS  test_vectorialbasisfunctionset)

add_executable(test_vectorialbasisfunctionset test-vectorialbasisfunctionset.cc)
dune_target_link_libraries(test_vectorialbasisfunctionset "${DUNE_LIBS};${LOCAL_LIBS}")
add_dune_mpi_flags(test_vectorialbasisfunctionset )

# include not needed for CMake
# include $(top_srcdir)/am/global-rulesset(TESTS  parametertest)
foreach(i ${TESTS})
  add_test(${i} ${i})
endforeach(i ${TESTS})
# We do not want want to build the tests during make all
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${TESTS} ${COMPILE_XFAIL_TESTS})
