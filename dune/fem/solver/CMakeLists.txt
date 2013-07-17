# don't follow the full GNU-standard; we definitely need automake >= 1.5
set(AUTOMAKE_OPTIONS  foreign 1.5)
set(SUBDIRS  ode oemsolver test)
set(noinst_LTLIBRARIES  libsolver.la)
# construct only from sub-libraries
add_executable(libsolver_la )
# if par DG was found, then ODE solver are not included
# if PARDG
# libsolver_la_LIBADD =
# else
set(libsolver_la_LIBADD  ode/libodesolver.la)
# endif
set(solverdir   ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/solver)
set(solver_HEADERS  cginverseoperator.hh
                 diagonalpreconditioner.hh
                 inverseoperators.hh
                 istlsolver.hh
                 odesolver.hh
                 odesolverinterface.hh
                 oemsolver.hh
                 multistep.hh
                 newtoninverseoperator.hh
                 pardg.hh
                 petscsolver.hh
                 rungekutta.hh
                 timeprovider.hh
                 umfpacksolver.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${solver_HEADERS} DESTINATION ${solverdir})
foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt EXCLUDE_FROM_ALL)
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
endforeach(i ${SUBDIRS})
