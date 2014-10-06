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

set(EXTRA_DIST  checkleafcodim1.hh)
set(LDADD  ${ALL_PKG_LDFLAGS} ${ALL_PKG_LIBS} ${LOCAL_LIBS})

set(AM_CPPFLAGS  ${ALL_PKG_CPPFLAGS})
set(TESTS  cachingquadrature)
set(check_PROGRAMS  cachingquadrature)
add_executable(cachingquadrature cachingquadrature.cc)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
foreach(i ${TESTS})
  add_test(${i} ${i})
endforeach(i ${TESTS})
# We do not want want to build the tests during make all
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${TESTS} ${COMPILE_XFAIL_TESTS})

