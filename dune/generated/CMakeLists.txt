add_library(generated_module SHARED EXCLUDE_FROM_ALL generated_module.cc)
target_include_directories(generated_module PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
set_target_properties(generated_module PROPERTIES PREFIX "")

configure_file(__init__.py __init__.py COPYONLY)

install(FILES generated_module.cc DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/python/dune/generated)
