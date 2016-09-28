set(modules
  AddPAPIFlags.cmake
  AddSIONlibFlags.cmake
  CommandLineHacks.cmake
  DuneFemMacros.cmake
  FemShort.cmake
  FindPAPI.cmake
  FindPThreads.cmake
  FindSIONlib.cmake
  FindXDR.cmake
)
install(FILES ${modules} DESTINATION ${DUNE_INSTALL_MODULEDIR})
