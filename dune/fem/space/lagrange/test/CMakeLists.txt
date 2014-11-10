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

# copy data to build source to make tests work         
configure_file(2dgrid.dgf ${CMAKE_CURRENT_BINARY_DIR}/2dgrid.dgf COPYONLY)
configure_file(3dgrid.dgf ${CMAKE_CURRENT_BINARY_DIR}/3dgrid.dgf COPYONLY)

add_executable(tester main.cc)
dune_target_link_libraries( tester "${DUNE_LIBS};${LOCAL_LIBS}")
add_executable(generic generic.cc)
dune_target_link_libraries( generic "${DUNE_LIBS};${LOCAL_LIBS}")
add_executable(testnonscalarmapper testnonscalarmapper.cc)
dune_target_link_libraries( testnonscalarmapper "${DUNE_LIBS};${LOCAL_LIBS}" )

dune_add_test(tester generic testnonscalarmapper)
