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
add_definitions( "-DGRIDDIM=${GRIDDIM}" ) 
add_definitions( "-DDIMRANGE=${DIMRANGE}" )
add_definitions( "-DPOLORDER=${POLORDER}" )
add_definitions( "-DWANT_GRAPE=${WANT_GRAPE}" )

set(TESTS  lagrangeadapt 
           adapt 
           padapt
           dgcomm_4.sh 
           dgcomm_cached_4.sh
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
add_executable(adapt adapt.cc)
dune_target_link_libraries(adapt "${DUNE_LIBS};${LOCAL_LIBS}")
add_executable(padapt padapt.cc)
dune_target_link_libraries(padapt "${DUNE_LIBS};${LOCAL_LIBS}")
add_executable(combinedspacetest combinedspacetest.cc)
dune_target_link_libraries(combinedspacetest "${DUNE_LIBS};${LOCAL_LIBS}")
add_executable(lagrangeinterpolation lagrangeinterpolation.cc)
dune_target_link_libraries(lagrangeinterpolation "${DUNE_LIBS};${LOCAL_LIBS}")
add_executable(interpolation interpolation.cc)
dune_target_link_libraries(interpolation "${DUNE_LIBS};${LOCAL_LIBS}")
add_executable(lagrangeadapt lagrangeadapt.cc)
dune_target_link_libraries(lagrangeadapt "${DUNE_LIBS};${LOCAL_LIBS}")
add_executable(lagrangeglobalrefine lagrangeglobalrefine.cc)
dune_target_link_libraries(lagrangeglobalrefine "${DUNE_LIBS};${LOCAL_LIBS}")
add_executable(legendretest legendretest.cc)
dune_target_link_libraries(legendretest "${DUNE_LIBS};${LOCAL_LIBS}")
add_executable(arraytest arraytest.cc)
dune_target_link_libraries(arraytest "${DUNE_LIBS};${LOCAL_LIBS}")
add_executable(dgcomm dgcomm.cc)
dune_target_link_libraries(dgcomm "${DUNE_LIBS};${LOCAL_LIBS}")
add_executable(dgcomm_cached dgcomm.cc)
dune_target_link_libraries(dgcomm_cached "${DUNE_LIBS};${LOCAL_LIBS}")

add_dune_alugrid_flags(lagrangeadapt)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules
foreach(i ${TESTS})
  add_test(${i} ${i})
endforeach(i ${TESTS})
# We do not want want to build the tests during make all
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${TESTS} ${COMPILE_XFAIL_TESTS})
