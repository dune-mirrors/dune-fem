set(modules 
  DuneFemMacros.cmake
  FindPAPI.cmake
  FindPThreads.cmake
  FindSIONlib.cmake
)
install(FILES ${modules} DESTINATION ${DUNE_INSTALL_MODULEDIR})
