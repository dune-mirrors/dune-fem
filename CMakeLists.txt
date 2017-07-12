add_subdirectory(dune)

configure_file(setup.py.in setup.py)
if(NOT (${DUNE_COMMON_VERSION} VERSION_LESS 2.6))
  configure_file(setup.py.in setup.py)
endif()
