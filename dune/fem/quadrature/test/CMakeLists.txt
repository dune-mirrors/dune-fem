# set default grid type for make all
set(GRIDTYPE YASPGRID)
#set(GRIDTYPE ALUGRID_SIMPLEX)
#set(GRIDTYPE ALUGRID_CUBE)
#set(GRIDTYPE UGGRID)
set(GRIDDIM 3)

add_definitions( "-D${GRIDTYPE}" )
add_definitions( "-DGRIDDIM=${GRIDDIM}" ) 

# copy data to build source to make tests work         
configure_file(1dgrid.dgf ${CMAKE_CURRENT_BINARY_DIR}/1dgrid.dgf COPYONLY)
configure_file(2dgrid.dgf ${CMAKE_CURRENT_BINARY_DIR}/2dgrid.dgf COPYONLY)
configure_file(3dgrid.dgf ${CMAKE_CURRENT_BINARY_DIR}/3dgrid.dgf COPYONLY)
configure_file(3dcubesimp.dgf ${CMAKE_CURRENT_BINARY_DIR}/3dcubesimp.dgf COPYONLY)
configure_file(parameter ${CMAKE_CURRENT_BINARY_DIR}/parameter COPYONLY)

add_executable(cachingquadrature cachingquadrature.cc)
dune_target_link_libraries(cachingquadrature "${DUNE_LIBS};${LOCAL_LIBS}")

dune_fem_add_test(cachingquadrature)

