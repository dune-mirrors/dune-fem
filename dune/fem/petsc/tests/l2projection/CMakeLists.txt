# IGNORE_HEADERS
#GRIDTYPE = SPGRID
set(GRIDTYPE  YASPGRID)
#GRIDTYPE = ALUGRID_SIMPLEX
set(GRIDDIM  2)
set(POLORDER  1)
set(EXTRA_DIST  parameter unitcube1.dgf unitcube2.dgf unitcube3.dgf base.hh massoperator.hh
                                                 mpiexec-petscmatrix.sh mpiexec.sh)
# linker flags
set(LDADD  ${ALL_PKG_LDFLAGS} ${ALL_PKG_LIBS} ${DUNEMPILIBS} ${LOCAL_LIBS} ${MPI_LDFLAGS})
# preprocessor flags
set(AM_CPPFLAGS  ${ALL_PKG_CPPFLAGS} ${DUNEMPICPPFLAGS})
set(CHECK_PROGS  l2projection l2projection_petscmatrix)
set(TESTS  ${CHECK_PROGS})
set(check_PROGRAMS  ${CHECK_PROGS})
set(CLEANFILES  *.tex data/*)
add_executable(petsc_l2projection l2projection.cc)
target_link_libraries(petsc_l2projection ${DUNE_LIBS})
add_executable(l2projection_petscmatrix l2projection.cc)
target_link_libraries(l2projection_petscmatrix ${DUNE_LIBS})
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
foreach(i ${TESTS})
  add_test(${i} ${i})
endforeach(i ${TESTS})
# We do not want want to build the tests during make all
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${TESTS} ${COMPILE_XFAIL_TESTS})