# IGNORE_HEADERS
# set default grid type for make all
set(GRIDTYPE SGRID)
set(GRIDDIM 2)
# test programs
set(check_PROGRAMS  test-filteredgridpart lagrangeinterpolation test-idgridpart test-geogridpart)
set(LDADD  ${ALL_PKG_LDFLAGS} ${ALL_PKG_LIBS} ${LOCAL_LIBS})
set(AM_CPPFLAGS  ${ALL_PKG_CPPFLAGS})
add_executable(gp_lagrangeinterpolation lagrangeinterpolation.cc)
target_link_libraries(gp_lagrangeinterpolation ${DUNE_LIBS})
add_executable(test_filteredgridpart test-filteredgridpart.cc)
target_link_libraries(test_filteredgridpart ${DUNE_LIBS})
add_executable(test_idgridpart test-idgridpart.cc)
add_executable(test_geogridpart test-geogridpart.cc)
set(EXTRA_DIST  1dgrid.dgf 2dgrid.dgf 3dgrid.dgf
             checkgeometry.hh  checkindexset.hh
             checkintersections.hh  checkseed.hh
             failure.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
