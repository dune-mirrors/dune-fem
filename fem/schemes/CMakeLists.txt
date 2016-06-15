#install headers
install(FILES
  dgelliptic.hh dgrhs.hh dirichletconstraints.hh estimator.hh interpolationscheme.hh solver.hh dgfemscheme.hh diffusionmodel.hh elliptic.hh femscheme.hh rhs.hh
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/schemes)
