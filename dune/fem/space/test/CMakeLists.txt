set(POLORDER  2)
set(DIMRANGE  3)
set(WANT_GRAPE  0)
# set default grid type for make all
set(GRIDTYPE YASPGRID)
#set( GRIDTYPE PARALLELGRID_ALUGRID_SIMPLEX )
#set( GRIDTYPE ALUGRID_SIMPLEX )
#set( GRIDTYPE ALUGRID_CONFORM )
#set( GRIDTYPE ALUGRID_CUBE ) 
#set( GRIDTYPE ALBERTAGRID )
set(GRIDDIM 2)

add_definitions( "-D${GRIDTYPE}" )
add_definitions( "-DDIMRANGE=${DIMRANGE}" )
add_definitions( "-DPOLORDER=${POLORDER}" )
add_definitions( "-DWANT_GRAPE=${WANT_GRAPE}" )

set(TESTS  lagrangeadapt 
           adapt 
           padapt
           arraytest
           combinedspacetest
           dgcomm
           dgcomm_cached
           interpolation
           l2projection
           lagrangeglobalrefine
           lagrangeinterpolation
           legendretest)

add_executable(l2projection l2projection.cc)
dune_target_link_libraries(l2projection "${DUNE_LIBS};${LOCAL_LIBS}")
set_property(TARGET l2projection APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=${GRIDDIM}" )
add_executable(adapt adapt.cc)
dune_target_link_libraries(adapt "${DUNE_LIBS};${LOCAL_LIBS}")
set_property(TARGET adapt APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=${GRIDDIM}" )
add_executable(padapt padapt.cc)
dune_target_link_libraries(padapt "${DUNE_LIBS};${LOCAL_LIBS}")
set_property(TARGET padapt APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=${GRIDDIM}" )
add_executable(combinedspacetest combinedspacetest.cc)
dune_target_link_libraries(combinedspacetest "${DUNE_LIBS};${LOCAL_LIBS}")
set_property(TARGET combinedspacetest APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=${GRIDDIM}" )
add_executable(lagrangeinterpolation lagrangeinterpolation.cc)
dune_target_link_libraries(lagrangeinterpolation "${DUNE_LIBS};${LOCAL_LIBS}")
set_property(TARGET lagrangeinterpolation APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=${GRIDDIM}" )
add_executable(interpolation interpolation.cc)
dune_target_link_libraries(interpolation "${DUNE_LIBS};${LOCAL_LIBS}")
set_property(TARGET interpolation APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=${GRIDDIM}" )
add_executable(lagrangeadapt lagrangeadapt.cc)
dune_target_link_libraries(lagrangeadapt "${DUNE_LIBS};${LOCAL_LIBS}")
set_property(TARGET lagrangeadapt APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=${GRIDDIM}" )
add_executable(lagrangeglobalrefine lagrangeglobalrefine.cc)
dune_target_link_libraries(lagrangeglobalrefine "${DUNE_LIBS};${LOCAL_LIBS}")
set_property(TARGET lagrangeglobalrefine APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=${GRIDDIM}" )
add_executable(legendretest legendretest.cc)
dune_target_link_libraries(legendretest "${DUNE_LIBS};${LOCAL_LIBS}")
set_property(TARGET legendretest APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=${GRIDDIM}" )
add_executable(arraytest arraytest.cc)
dune_target_link_libraries(arraytest "${DUNE_LIBS};${LOCAL_LIBS}")
set_property(TARGET arraytest APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=${GRIDDIM}" )
add_executable(dgcomm dgcomm.cc)
dune_target_link_libraries(dgcomm "${DUNE_LIBS};${LOCAL_LIBS}")
set_property(TARGET dgcomm APPEND PROPERTY COMPILE_DEFINITIONS "UGRIDDIM" )
set_property(TARGET dgcomm APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=3" )
set_property(TARGET dgcomm APPEND PROPERTY COMPILE_DEFINITIONS "UWORLDDIM" )
set_property(TARGET dgcomm APPEND PROPERTY COMPILE_DEFINITIONS "WORLDDIM=3" )
set_property(TARGET dgcomm APPEND PROPERTY COMPILE_DEFINITIONS "WANT_CACHED_COMM_MANAGER=0" )
add_executable(dgcomm_cached dgcomm.cc)
dune_target_link_libraries(dgcomm_cached "${DUNE_LIBS};${LOCAL_LIBS}")
set_property(TARGET dgcomm_cached APPEND PROPERTY COMPILE_DEFINITIONS "UGRIDDIM" )
set_property(TARGET dgcomm_cached APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=3" )
set_property(TARGET dgcomm_cached APPEND PROPERTY COMPILE_DEFINITIONS "UWORLDDIM" )
set_property(TARGET dgcomm_cached APPEND PROPERTY COMPILE_DEFINITIONS "WORLDDIM=3" )



# todo: dgcomm_cached_4.sh and dgcomm_4.sh are not used to verify mpi tests
#add_test(dgcomm_4.sh dgcomm-4.sh)
#add_test(dgcomm_cached_4.sh dgcomm-cached-4.sh)

# We do not want want to build the tests during make all
# but just build them on demand

# include not needed for CMake
# include $(top_srcdir)/am/global-rules
foreach(i ${TESTS})
  add_test(${i} ${i})
endforeach(i ${TESTS})
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${TESTS} ${COMPILE_XFAIL_TESTS})
