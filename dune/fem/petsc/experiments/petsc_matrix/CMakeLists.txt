set(petsc_matrixdir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/petsc/experiments/petsc_matrix)
set(petsc_matrix_HEADERS  function.hh l2error.hh laplaceassembler.hh printer.hh
                                                                                         scheme.hh solver.hh)
#EXTRA_DIST = parameter unitcube1.dgf unitcube2.dgf unitcube3.dgf
# linker flags
set(LDADD  ${ALL_PKG_LDFLAGS} ${ALL_PKG_LIBS} ${DUNEMPILIBS} ${LOCAL_LIBS})
# preprocessor flags
set(AM_CPPFLAGS  ${ALL_PKG_CPPFLAGS} ${DUNEMPICPPFLAGS})
set(CHECK_PROGS  )
set(check_PROGRAMS  ${CHECK_PROGS})
set(CLEANFILES  *.tex data/*)

# add_executable(matrix_test matrix_test.cc)
# dune_target_link_libraries(matrix_test ${DUNE_LIBS})

# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${petsc_matrix_HEADERS} DESTINATION ${petsc_matrixdir})