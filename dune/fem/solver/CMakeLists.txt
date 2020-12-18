dune_install(cginverseoperator.hh diagonalpreconditioner.hh
             inverseoperatorinterface.hh istl.hh istlinverseoperators.hh
             krylovinverseoperators.hh multistep.hh newtoninverseoperator.hh
             odesolver.hh odesolverinterface.hh oemsolver.hh
             parameter.hh pardginverseoperators.hh
             petscinverseoperators.hh preconditionedinverseoperator.hh
             timeprovider.hh umfpacksolver.hh spqrsolver.hh ldlsolver.hh eigen.hh)

dune_add_subdirs(communication linear rungekutta test)

dune_add_library(solver OBJECT
    rungekutta/butchertable.cc )

# avoid deprecation warning.
exclude_from_headercheck(
  oemsolver.hh
  pardginverseoperators.hh
  )
