add_subdirectory(common)
add_subdirectory(filter)
add_subdirectory(geometrygridpart)

install(FILES geometrygridpart.hh DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/gridpart)
