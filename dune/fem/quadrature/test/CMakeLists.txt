# IGNORE_HEADERS
#CXXFLAGS=-g -Wall
# set default grid type for make all
set(GRIDTYPE YASPGRID)
#set(GRIDTYPE PRISMGRID_ALUGRID_SIMPLEX)
#set(GRIDTYPE ALUGRID_SIMPLEX)
#set(GRIDTYPE ALUGRID_CUBE)
#set(GRIDTYPE UGGRID)
set(GRIDDIM 3)

add_definitions( "-D${GRIDTYPE}" )
add_definitions( "-DGRIDDIM=${GRIDDIM}" ) 

set(TESTS  cachingquadrature)

# copy data to build source to make tests work         
configure_file(1dgrid.dgf ${CMAKE_CURRENT_BINARY_DIR}/1dgrid.dgf COPYONLY)
configure_file(2dgrid.dgf ${CMAKE_CURRENT_BINARY_DIR}/2dgrid.dgf COPYONLY)
configure_file(3dgrid.dgf ${CMAKE_CURRENT_BINARY_DIR}/3dgrid.dgf COPYONLY)
configure_file(3dcubesimp.dgf ${CMAKE_CURRENT_BINARY_DIR}/3dcubesimp.dgf COPYONLY)
configure_file(parameter ${CMAKE_CURRENT_BINARY_DIR}/parameter COPYONLY)

add_executable(cachingquadrature cachingquadrature.cc)
dune_target_link_libraries( cachingquadrature "${DUNE_LIBS};${LOCAL_LIBS}")

# include not needed for CMake
# include $(top_srcdir)/am/global-rules
foreach(i ${TESTS})
  add_test(${i} ${i})
endforeach(i ${TESTS})
# We do not want want to build the tests during make all
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${TESTS} ${COMPILE_XFAIL_TESTS})

