# IGNORE_HEADERS
# set default grid type for make all
set(GRIDTYPE SGRID)
set(GRIDDIM 2)

add_definitions( "-D${GRIDTYPE}" )
add_definitions( "-DGRIDDIM=${GRIDDIM}" ) 
# copy data to build source to make tests work         


configure_file(1dgrid.dgf ${CMAKE_CURRENT_BINARY_DIR}/1dgrid.dgf COPYONLY)
configure_file(2dgrid.dgf ${CMAKE_CURRENT_BINARY_DIR}/2dgrid.dgf COPYONLY)
configure_file(3dgrid.dgf ${CMAKE_CURRENT_BINARY_DIR}/3dgrid.dgf COPYONLY)

# test programs
set(TESTS  gp_lagrangeinterpolation 
           test_idgridpart 
           test_geogridpart)

# non compiling tests
# todo: fix compiler errors
#set(TESTS "${TESTS} test_filteredgridpart")


add_executable(gp_lagrangeinterpolation lagrangeinterpolation.cc)
dune_target_link_libraries(gp_lagrangeinterpolation "${DUNE_LIBS};${LOCAL_LIBS}")
add_executable(test_filteredgridpart test-filteredgridpart.cc)
dune_target_link_libraries(test_filteredgridpart "${DUNE_LIBS};${LOCAL_LIBS}")
add_executable(test_idgridpart test-idgridpart.cc)
dune_target_link_libraries(test_idgridpart "${DUNE_LIBS};${LOCAL_LIBS}")
add_executable(test_geogridpart test-geogridpart.cc)
dune_target_link_libraries(test_geogridpart "${DUNE_LIBS};${LOCAL_LIBS}")

set_property(TARGET gp_lagrangeinterpolation APPEND PROPERTY 
             COMPILE_DEFINITIONS "COUNT_FLOPS;USE_FILTEREDGRIDPART;POLORDER=1" )
set_property(TARGET test_filteredgridpart APPEND PROPERTY 
             COMPILE_DEFINITIONS "COUNT_FLOPS;USE_FILTEREDGRIDPART" )


# include not needed for CMake
# include $(top_srcdir)/am/global-rules
foreach(i ${TESTS})
  add_test(${i} ${i})
endforeach(i ${TESTS})
# We do not want want to build the tests during make all
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${TESTS} ${COMPILE_XFAIL_TESTS})

