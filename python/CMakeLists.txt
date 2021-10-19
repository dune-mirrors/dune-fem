add_subdirectory(dune)

dune_python_install_package(
  PATH "."
  CMAKE_METADATA_FILE dune/fem/metadata.cmake
  DEPENDS _fem
  CMAKE_METADATA_FLAGS DUNE_OPTS_FILE
)
if ("${CMAKE_VERSION}" VERSION_GREATER 3.14)
  cmake_policy(SET CMP0087 NEW)
endif ("${CMAKE_VERSION}" VERSION_GREATER 3.14)
