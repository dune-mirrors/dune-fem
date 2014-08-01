set(commtestdir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/petsc/experiments/commtest)
#GRIDTYPE = ALUGRID_SIMPLEX
#GRIDTYPE = ALUGRID_CUBE
#GRIDTYPE = SPGRID
set(GRIDTYPE  YASPGRID)
#YASPGRID = 1
set(GRIDDIM  2)
set(DIMRANGE  2)
set(POLORDER  1)
set(AM_CPPFLAGS  ${ALL_PKG_CPPFLAGS} -DFEMTIMER -DGRIDDIM=${GRIDDIM} -DDIMRANGE=${DIMRANGE}
                          -DPOLORDER=${POLORDER})
set(LDADD  ${ALL_PKG_LDFLAGS} ${ALL_PKG_LIBS} ${LOCAL_LIBS} ${MPI_LDFLAGS})
set(CHECK_PROGS  )
set(check_PROGRAMS  ${CHECK_PROGS})
set(CLEANFILES  *.tex data/*)

# add_executable(commtest commtest.cc)
# dune_target_link_libraries(commtest ${DUNE_LIBS})

set(EXTRA_DIST  domain_commtest.dgf parameter_commtest)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
