set(SUBDIRS  ode oemsolver test)

# construct only from sub-libraries
# if par DG was found, then ODE solver are not included
# if PARDG
# libsolver_la_LIBADD =
# else

# endif
set(solverdir   ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/solver)
set(solver_HEADERS cginverseoperator.hh
                   diagonalpreconditioner.hh
                   istlinverseoperators.hh
                   istlsolver.hh
                   multistep.hh
                   newtoninverseoperator.hh
                   odesolver.hh
                   odesolverinterface.hh
                   oemsolver.hh
                   pardg.hh
                   pardginverseoperators.hh
                   petscsolver.hh
                   preconditionedinverseoperator.hh
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

dune_add_library(solver OBJECT 
    ode/odeobjectfiles.cc 
    ode/odeobjectfiles_mpi.cc
  )
