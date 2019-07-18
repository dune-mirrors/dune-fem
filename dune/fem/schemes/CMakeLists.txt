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

if( NOT dune-python_FOUND )
  exclude_from_headercheck( ${HEADERS} )
endif()

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/schemes)
