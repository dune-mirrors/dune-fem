set(HEADERS
  conservationlawmodel.hh
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
  molgalerkin.hh
)

exclude_from_headercheck( dgfemscheme.hh )

if( NOT DUNE_ENABLE_PYTHONBINDINGS )
  exclude_from_headercheck( ${HEADERS} )
endif()

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/schemes)
