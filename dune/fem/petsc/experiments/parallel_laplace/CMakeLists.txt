set(parallel_laplacedir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/petsc/experiments/parallel_laplace)
set(parallel_laplace_HEADERS  function.hh l2error.hh laplaceassembler.hh printer.hh
                                                                                                  scheme.hh solver.hh)
# linker flags
set(LDADD  ${ALL_PKG_LDFLAGS} ${ALL_PKG_LIBS} ${DUNEMPILIBS} ${LOCAL_LIBS})
# preprocessor flags
set(AM_CPPFLAGS  ${ALL_PKG_CPPFLAGS} ${DUNEMPICPPFLAGS})
set(CHECK_PROGS  )
set(check_PROGRAMS  ${CHECK_PROGS})
set(CLEANFILES  *.tex data/* rhs.out.gnuplot* solution.out.gnuplot* error.out.gnuplot*)
add_executable(parallel_laplace parallel_laplace.cc)
target_link_libraries(parallel_laplace ${DUNE_LIBS})
set(EXTRA_DIST  mpiexec.sh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${parallel_laplace_HEADERS} DESTINATION ${parallel_laplacedir})