set(LDADD  ${ALL_PKG_LDFLAGS} ${ALL_PKG_LIBS} ${LOCAL_LIBS})
set(AM_CPPFLAGS  ${ALL_PKG_CPPFLAGS})
# set default grid type for make all
set(GRIDTYPE YASPGRID)
set(GRIDDIM 2)
set(DIMENSION  2)
set(POLORDER  2)
set(GEOMETRYTYPE  simplex)
set(SCALAR  0)

add_definitions( "-D${GRIDTYPE}" )
add_definitions( "-DGRIDDIM=${GRIDDIM}" ) 
add_definitions( "-DDIMENSION=${DIMENSION}" )
add_definitions( "-DPOLORDER=${POLORDER}" )
add_definitions( "-DGEOMETRYTYPE=${GEOMETRYTYPE}" )
add_definitions( "-DSCALAR=${SCALAR}" )

set(TESTS  tester generic )
#testnonscalarmapper)
set(check_PROGRAMS  tester generic testnonscalarmapper)

add_executable(tester main.cc)
add_executable(generic generic.cc)
add_executable(testnonscalarmapper testnonscalarmapper.cc)
set(testdir ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/space/lagrange/test)
set(test_HEADERS  basisfunctiontest.hh  mappertest.hh basisfunctiontest.cc mappertest.cc)
set(EXTRA_DIST  2dgrid.dgf cube.dgf onesimplex.dgf)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules
foreach(i ${TESTS})
  add_test(${i} ${i})
endforeach(i ${TESTS})
# We do not want want to build the tests during make all
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${TESTS} ${COMPILE_XFAIL_TESTS})
