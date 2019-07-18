add_subdirectory(common)
add_subdirectory(function)
add_subdirectory(geometry)
add_subdirectory(grid)
add_subdirectory(space)
add_subdirectory(py)
add_subdirectory(pybind11)
add_subdirectory(quadrature)

set(HEADERS
  parameter.hh
)

if( NOT dune-python_FOUND )
  exclude_from_headercheck( ${HEADERS} )
endif()


install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/fempy)
