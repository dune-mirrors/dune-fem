add_subdirectory(dune)

dune_python_configure_bindings(
  PATH "."
  CMAKE_METADATA_FLAGS DUNE_OPTS_FILE PETSC_DIR
  REQUIRED "${DUNE_PYTHONBINDINGS_REQUIRED}"
)
if(POLICY CMP0087)
  cmake_policy(SET CMP0087 NEW)
endif()
