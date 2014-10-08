set( POLORD  3 )
set( GRIDTYPE  YASPGRID )
set( GRIDDIM  2 )

set( LDADD  ${ALL_PKG_LDFLAGS} ${ALL_PKG_LIBS} ${DUNEMPILIBS} ${LOCAL_LIBS})
set(AM_CPPFLAGS  ${ALL_PKG_CPPFLAGS} ${DUNEMPICPPFLAGS} -DPOLORDER=${POLORD})

# remove this until this test is fixed again
#set( ADAPTPROG  prod )
# remove if-condition after test is fixed again
if( ${ADAPTPROG} )
  add_definitions( "-D${GRIDTYPE}" )
  add_definitions( "-DGRIDDIM=${GRIDDIM}" ) 
  add_definitions( "-DPOLORD=${POLORD}" )
  
  set(TESTS  ${ADAPTPROG})
  set(check_PROGRAMS  ${ADAPTPROG})
  add_executable(prod productfunctiontest.cc)
  dune_target_link_libraries(prod ${DUNE_LIBS} ${LOCAL_LIBS})
  
  # include not needed for CMake
  # include $(top_srcdir)/am/global-rules
  foreach(i ${TESTS})
    add_test(${i} ${i})
  endforeach(i ${TESTS})
  # We do not want want to build the tests during make all
  # but just build them on demand
  add_directory_test_target(_test_target)
  add_dependencies(${_test_target} ${TESTS} ${COMPILE_XFAIL_TESTS})
endif( ${ADAPTPROG} )
