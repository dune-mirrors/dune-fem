# IGNORE_HEADERS
# set default grid type for make all
set(GRIDTYPE YASPGRID)
set(GRIDDIM 2)

add_definitions( "-D${GRIDTYPE}" )
add_definitions( "-DGRIDDIM=${GRIDDIM}" ) 

# also distribute the headers for the checks
set(TESTS  quad_tester)

add_executable(quad_tester main.cc refelem_test.cc cache_test.cc twist_test.cc point_test.cc)
dune_target_link_libraries(quad_tester "${DUNE_LIBS};${LOCAL_LIBS}")

dune_add_test(${TESTS})
