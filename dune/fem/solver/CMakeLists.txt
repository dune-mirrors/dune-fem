dune_install(cginverseoperator.hh diagonalpreconditioner.hh fempreconditioning.hh
             inverseoperatorinterface.hh istl.hh istlinverseoperators.hh
             krylovinverseoperators.hh multistep.hh newtoninverseoperator.hh
             odesolver.hh odesolverinterface.hh oemsolver.hh
             parameter.hh parameterdoc.hh pardginverseoperators.hh preconditionfunctionwrapper.hh
             petscinverseoperators.hh preconditionedinverseoperator.hh
             timeprovider.hh umfpacksolver.hh spqrsolver.hh ldlsolver.hh eigen.hh)

dune_add_subdirs(communication linear rungekutta test)

target_sources(dunefem PRIVATE
    rungekutta/butchertable.cc )

# avoid deprecation warning.
exclude_from_headercheck(
  oemsolver.hh
  pardginverseoperators.hh
  )
