set(HEADERS
  dgelliptic.hh
  dgfemscheme.hh
  dgrhs.hh
  diffusionmodel.hh
  dirichletconstraints.hh
  dirichletwrapper.hh
  elliptic.hh
  estimator.hh
  femscheme.hh
  galerkin.hh
  integrands.hh
  linearized.hh
  rhs.hh
)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/schemes)
