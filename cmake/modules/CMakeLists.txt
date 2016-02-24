set(modules 
  DuneFemMacros.cmake
  FindPAPI.cmake
  FindPThreads.cmake
  FindSIONlib.cmake
)
install(FILES "${modules}" DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/cmake/modules)
