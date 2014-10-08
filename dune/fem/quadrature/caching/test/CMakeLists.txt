# IGNORE_HEADERS
# set default grid type for make all
set(GRIDTYPE YASPGRID)
set(GRIDDIM 2)

add_definitions( "-D${GRIDTYPE}" )
add_definitions( "-DGRIDDIM=${GRIDDIM}" ) 

# also distribute the headers for the checks
set(TESTS  tester)

add_executable(quad_tester main.cc ../../../misc/test.cc ../../../misc/suite.cc refelem_test.cc cache_test.cc twist_test.cc point_test.cc)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules
foreach(i ${TESTS})
  add_test(${i} ${i})
endforeach(i ${TESTS})
# We do not want want to build the tests during make all
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${TESTS} ${COMPILE_XFAIL_TESTS})
