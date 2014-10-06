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

set(LDADD  ${ALL_PKG_LDFLAGS} ${ALL_PKG_LIBS} ${LOCAL_LIBS} ${MPI_LDFLAGS})
set(AM_CPPFLAGS  ${ALL_PKG_CPPFLAGS} -DPOLORDER=${POLORDER} -DDIMRANGE=${DIMRANGE} -DWANT_GRAPE=${WANT_GRAPE})
set(ADAPTPROG  lagrangeadapt adapt padapt)
# programs run during make check
set(TESTS  ${ADAPTPROG} dgcomm-4.sh dgcomm-cached-4.sh)
set(check_PROGRAMS  ${ADAPTPROG}
                 arraytest
                 combinedspacetest
                 dgcomm
                 dgcomm_cached
                 interpolation
     l2projection
     lagrangeglobalrefine
                 lagrangeinterpolation
                 legendretest)
set(EXTRA_DIST  systemmatrix.cc)

add_executable(l2projection l2projection.cc)
add_executable(adapt adapt.cc)
add_executable(padapt padapt.cc)
add_executable(combinedspacetest combinedspacetest.cc)
add_executable(lagrangeinterpolation lagrangeinterpolation.cc)
add_executable(interpolation interpolation.cc)
add_executable(lagrangeadapt lagrangeadapt.cc)
add_executable(lagrangeglobalrefine lagrangeglobalrefine.cc)
add_executable(legendretest legendretest.cc)
add_executable(arraytest arraytest.cc)
add_executable(dgcomm dgcomm.cc)
add_executable(dgcomm_cached dgcomm.cc)

dune_target_link_libraries(dgcomm ${DUNE_LIBS})
dune_target_link_libraries(dgcomm_cached ${DUNE_LIBS})

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
