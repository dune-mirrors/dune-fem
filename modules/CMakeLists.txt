SET(MODULES
  DuneFempyMacros.cmake
  UseModelCompiler.cmake
)

install(FILES ${MODULES} DESTINATION ${DUNE_INSTALL_MODULEDIR})
