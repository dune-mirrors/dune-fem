set(POLORDER  2)
set(DIMRANGE  3)
set(WANT_GRAPE  0)
#set( GRIDTYPE PARALLELGRID_ALUGRID_SIMPLEX )
#set( GRIDTYPE ALUGRID_SIMPLEX )
#set( GRIDTYPE ALUGRID_CONFORM )
#set( GRIDTYPE ALUGRID_CUBE )
#set( GRIDTYPE ALBERTAGRID )
set(GRIDDIM 2)
set(GRIDTYPE YASPGRID)

add_definitions( "-D${GRIDTYPE}" )
add_definitions( "-DDIMRANGE=${DIMRANGE}" )
add_definitions( "-DPOLORDER=${POLORDER}" )
add_definitions( "-DWANT_GRAPE=${WANT_GRAPE}" )

# copy data to build source to make tests work
configure_file(2dgrid.dgf ${CMAKE_CURRENT_BINARY_DIR}/2dgrid.dgf COPYONLY)
configure_file(3dgrid.dgf ${CMAKE_CURRENT_BINARY_DIR}/3dgrid.dgf COPYONLY)
configure_file(parameter ${CMAKE_CURRENT_BINARY_DIR}/parameter COPYONLY)
configure_file(dgcomm-4.sh ${CMAKE_CURRENT_BINARY_DIR}/dgcomm_4.sh COPYONLY)
configure_file(dgcomm-cached-4.sh ${CMAKE_CURRENT_BINARY_DIR}/dgcomm_cached_4.sh COPYONLY)

add_executable(l2projection l2projection.cc)
dune_target_link_libraries(l2projection "${DUNE_LIBS};${LOCAL_LIBS}")
set_property(TARGET l2projection APPEND PROPERTY
             COMPILE_DEFINITIONS "GRIDDIM=${GRIDDIM}" )
add_executable(adapt adapt.cc)
dune_target_link_libraries(adapt "${DUNE_LIBS};${LOCAL_LIBS}")
set_property(TARGET adapt APPEND PROPERTY
             COMPILE_DEFINITIONS "GRIDDIM=${GRIDDIM}" )
add_executable(padapt padapt.cc)
dune_target_link_libraries(padapt "${DUNE_LIBS};${LOCAL_LIBS}")
set_property(TARGET padapt APPEND PROPERTY
             COMPILE_DEFINITIONS "GRIDDIM=${GRIDDIM}" )
add_executable(combinedspacetest combinedspacetest.cc)
dune_target_link_libraries(combinedspacetest "${DUNE_LIBS};${LOCAL_LIBS}")
set_property(TARGET combinedspacetest APPEND PROPERTY
             COMPILE_DEFINITIONS "GRIDDIM=${GRIDDIM}" )
add_executable(lagrangeinterpolation lagrangeinterpolation.cc)
dune_target_link_libraries(lagrangeinterpolation "${DUNE_LIBS};${LOCAL_LIBS}")
set_property(TARGET lagrangeinterpolation APPEND PROPERTY
             COMPILE_DEFINITIONS "GRIDDIM=${GRIDDIM}" )
add_executable(interpolation interpolation.cc)
dune_target_link_libraries(interpolation "${DUNE_LIBS};${LOCAL_LIBS}")
set_property(TARGET interpolation APPEND PROPERTY
             COMPILE_DEFINITIONS "GRIDDIM=${GRIDDIM}" )
add_executable(hierarchicspace hierarchicspace.cc)
dune_target_link_libraries(hierarchicspace "${DUNE_LIBS};${LOCAL_LIBS}")
set_property(TARGET hierarchicspace APPEND PROPERTY
             COMPILE_DEFINITIONS "GRIDDIM=${GRIDDIM}" )
add_executable(lagrangeadapt lagrangeadapt.cc)
dune_target_link_libraries(lagrangeadapt "${DUNE_LIBS};${LOCAL_LIBS}")
set_property(TARGET lagrangeadapt APPEND PROPERTY
             COMPILE_DEFINITIONS "GRIDDIM=${GRIDDIM}" )
add_executable(lagrangeglobalrefine lagrangeglobalrefine.cc)
dune_target_link_libraries(lagrangeglobalrefine "${DUNE_LIBS};${LOCAL_LIBS}")
set_property(TARGET lagrangeglobalrefine APPEND PROPERTY
             COMPILE_DEFINITIONS "GRIDDIM=${GRIDDIM}" )
add_executable(legendretest legendretest.cc)
dune_target_link_libraries(legendretest "${DUNE_LIBS};${LOCAL_LIBS}")
set_property(TARGET legendretest APPEND PROPERTY
             COMPILE_DEFINITIONS "GRIDDIM=${GRIDDIM}" )
add_executable(localadapter localadapter.cc)
dune_target_link_libraries(localadapter "${DUNE_LIBS};${LOCAL_LIBS}")
set_property(TARGET localadapter APPEND PROPERTY
             COMPILE_DEFINITIONS "GRIDDIM=${GRIDDIM}" )
add_executable(arraytest arraytest.cc)
dune_target_link_libraries(arraytest "${DUNE_LIBS};${LOCAL_LIBS}")
set_property(TARGET arraytest APPEND PROPERTY
             COMPILE_DEFINITIONS "GRIDDIM=${GRIDDIM}" )
add_executable(dgcomm dgcomm.cc)
dune_target_link_libraries(dgcomm "${DUNE_LIBS};${LOCAL_LIBS}")
set_property(TARGET dgcomm APPEND PROPERTY
             COMPILE_DEFINITIONS "GRIDDIM=3;WORLDDIM=3;WANT_CACHED_COMM_MANAGER=0" )
add_executable(dgcomm_cached dgcomm.cc)
dune_target_link_libraries(dgcomm_cached "${DUNE_LIBS};${LOCAL_LIBS}")
set_property(TARGET dgcomm_cached APPEND PROPERTY
             COMPILE_DEFINITIONS "GRIDDIM=3;WORLDDIM=3" )

dune_add_test(lagrangeadapt adapt padapt
              arraytest combinedspacetest interpolation hierarchicspace
              l2projection lagrangeglobalrefine lagrangeinterpolation
              legendretest localadapter
              NO_DEPENDENCY dgcomm_4.sh dgcomm_cached_4.sh
              DEPENDENCY_ONLY dgcomm dgcomm_cached )
