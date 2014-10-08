set(POLORDER  2)
set(WANT_GRAPE  0)
set(GRIDTYPE YASPGRID)
set(GRIDDIM 2)

add_definitions( "-D${GRIDTYPE}" )
add_definitions( "-DGRIDDIM=${GRIDDIM}" ) 
add_definitions( "-DPOLORDER=${POLORDER}" )
add_definitions( "-DWANT_GRAPE=${WANT_GRAPE}" )


set(LDADD  ${ALL_PKG_LDFLAGS} ${ALL_PKG_LIBS} ${LOCAL_LIBS} ${MPI_LDFLAGS})
set(check_PROGRAMS  persistencemanagertest dataoutputtest)
set(TESTS  persistencemanagertest dataoutputtest)

add_executable(persistencemanagertest persistencemanagertest.cc)
dune_target_link_libraries(persistencemanagertest ${DUNE_LIBS} ${LOCAL_LIBS})
add_executable(dataoutputtest dataoutputtest.cc)
dune_target_link_libraries(dataoutputtest ${DUNE_LIBS} ${LOCAL_LIBS})

set(AM_CPPFLAGS  ${ALL_PKG_CPPFLAGS} -DPOLORDER=${POLORDER} -DWANT_GRAPE=${WANT_GRAPE})
set(CXXFLAGS  -g)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules
foreach(i ${TESTS})
  add_test(${i} ${i})
endforeach(i ${TESTS})
# We do not want want to build the tests during make all
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${TESTS} ${COMPILE_XFAIL_TESTS})
