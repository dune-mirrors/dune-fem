# linker flags
set(LDADD  ${ALL_PKG_LDFLAGS} ${ALL_PKG_LIBS} ${DUNEMPILIBS} ${LOCAL_LIBS})
# preprocessor flags
set(AM_CPPFLAGS  ${ALL_PKG_CPPFLAGS} ${DUNEMPICPPFLAGS})
set(EXTRA_DIST  testparameterfile
                                                 oldparameterfile
                                                 param_DarcyStokes
                                                 param_NavierStokes
                                                 date.sh
                                                 bc.sh)
set(TESTS  parametertest)
set(check_PROGRAMS  parametertest)
add_executable(parametertest parametertest.cc)
dune_target_link_libraries(parametertest ${DUNE_LIBS} ${LOCAL_LIBS})
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
foreach(i ${TESTS})
  add_test(${i} ${i})
endforeach(i ${TESTS})
# We do not want want to build the tests during make all
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${TESTS} ${COMPILE_XFAIL_TESTS})
