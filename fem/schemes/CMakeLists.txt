set(HEADERS
  dgelliptic.hh
  dgfemscheme.hh
  dgrhs.hh
  diffusionmodel.hh
  dirichletconstraints.hh
  elliptic.hh
  estimator.hh
  femscheme.hh
  galerkin.hh
  rhs.hh
  solver.hh
)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/schemes)
