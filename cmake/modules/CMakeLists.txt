set(modules
  AddPAPIFlags.cmake
  AddSIONlibFlags.cmake
  CommandLineHacks.cmake
  CorrectWindowsPaths.cmake
  DuneFemMacros.cmake
  FemShort.cmake
  FindPackageMultipass.cmake
  FindPAPI.cmake
  FindPETSc.cmake
  FindPThreads.cmake
  FindSIONlib.cmake
  FindXDR.cmake
  ResolveCompilerPaths.cmake
)
install(FILES ${modules} DESTINATION ${DUNE_INSTALL_MODULEDIR})
