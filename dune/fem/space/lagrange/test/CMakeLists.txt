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

set(TESTS  tester 
           generic 
           testnonscalarmapper)

add_executable(tester main.cc)
dune_target_link_libraries( tester ${DUNE_LIBS} ${LOCAL_LIBS})
add_executable(generic generic.cc)
dune_target_link_libraries( generic ${DUNE_LIBS}  ${LOCAL_LIBS})
add_executable(testnonscalarmapper testnonscalarmapper.cc)
dune_target_link_libraries( testnonscalarmapper ${DUNE_LIBS} ${LOCAL_LIBS} )

# include not needed for CMake
# include $(top_srcdir)/am/global-rules
foreach(i ${TESTS})
  add_test(${i} ${i})
endforeach(i ${TESTS})
# We do not want want to build the tests during make all
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${TESTS} ${COMPILE_XFAIL_TESTS})
