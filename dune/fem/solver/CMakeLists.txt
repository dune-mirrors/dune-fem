dune_install(cginverseoperator.hh diagonalpreconditioner.hh istlinverseoperators.hh
             istlsolver.hh multistep.hh newtoninverseoperator.hh
             odesolver.hh odesolverinterface.hh oemsolver.hh
             pardg.hh pardginverseoperators.hh petscsolver.hh
             preconditionedinverseoperator.hh rungekutta.hh timeprovider.hh
             umfpacksolver.hh spqrsolver.hh)

dune_add_subdirs(rungekutta ode oemsolver test)

dune_add_library(solver OBJECT
    ode/odeobjectfiles.cc
    ode/odeobjectfiles_mpi.cc
    rungekutta/butchertable.cc )
