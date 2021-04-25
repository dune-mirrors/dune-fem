add_subdirectory(dune)

dune_python_install_package(
  PATH "."
  CMAKE_METADATA_FILE dune/fem/metadata.cmake
  DEPENDS _fem
)
