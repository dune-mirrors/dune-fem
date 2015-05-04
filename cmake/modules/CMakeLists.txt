set(modules 
  DuneFemMacros.cmake
  GridSelector.cmake
)
install(FILES "${modules}" DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/cmake/modules)
