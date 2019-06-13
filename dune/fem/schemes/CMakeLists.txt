set(HEADERS
  dgelliptic.hh
  dgfemscheme.hh
  diffusionmodel.hh
  dirichletconstraints.hh
  dirichletwrapper.hh
  elliptic.hh
  femscheme.hh
  galerkin.hh
  integrands.hh
  linearized.hh
)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/schemes)
