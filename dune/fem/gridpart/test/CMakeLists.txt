# IGNORE_HEADERS
# set default grid type for make all
set(GRIDTYPE SGRID)
set(GRIDDIM 2)

add_definitions( "-D${GRIDTYPE}" )
add_definitions( "-DGRIDDIM=${GRIDDIM}" ) 

# test programs
set(TESTS  test-filteredgridpart lagrangeinterpolation test-idgridpart test-geogridpart)
set(check_PROGRAMS  test-filteredgridpart lagrangeinterpolation test-idgridpart test-geogridpart)
set(LDADD  ${ALL_PKG_LDFLAGS} ${ALL_PKG_LIBS} ${LOCAL_LIBS})
set(AM_CPPFLAGS  ${ALL_PKG_CPPFLAGS})

add_executable(gp_lagrangeinterpolation lagrangeinterpolation.cc)
dune_target_link_libraries(gp_lagrangeinterpolation ${DUNE_LIBS} ${LOCAL_LIBS})
add_executable(test_filteredgridpart test-filteredgridpart.cc)
dune_target_link_libraries(test_filteredgridpart ${DUNE_LIBS} ${LOCAL_LIBS})
add_executable(test_idgridpart test-idgridpart.cc)
dune_target_link_libraries(test_idgridpart ${DUNE_LIBS} ${LOCAL_LIBS})
add_executable(test_geogridpart test-geogridpart.cc)
dune_target_link_libraries(test_geogridpart ${DUNE_LIBS} ${LOCAL_LIBS})
set(EXTRA_DIST  1dgrid.dgf 2dgrid.dgf 3dgrid.dgf
             checkgeometry.hh  checkindexset.hh
             checkintersections.hh  checkseed.hh
             failure.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
foreach(i ${TESTS})
  add_test(${i} ${i})
endforeach(i ${TESTS})
# We do not want want to build the tests during make all
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${TESTS} ${COMPILE_XFAIL_TESTS})

