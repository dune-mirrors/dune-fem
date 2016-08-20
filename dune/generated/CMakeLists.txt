add_definitions( "-DALUGRID_CONFORM=1" )
add_definitions( "-DPOLORDER=1" )
add_definitions( "-DGRIDDIM=2" )
add_definitions( "-DWORLDDIM=2" )
add_definitions( "-DWANT_ISTL=1" )
add_definitions( "-DEQUAL_ORDER=0" )
add_definitions( "-DFULLSCHEME=0" )
add_definitions( "-DMINIELEMENT=0" )
add_definitions( "-DENABLE_UMFPACK=0" )

add_library(generated_module SHARED EXCLUDE_FROM_ALL generated_module.cc)
target_include_directories(generated_module PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
set_target_properties(generated_module PROPERTIES PREFIX "")

configure_file(__init__.py __init__.py COPYONLY)

install(FILES generated_module.cc DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/python/dune/generated)

add_library(modelimpl SHARED EXCLUDE_FROM_ALL modelimpl.cc)
target_include_directories(modelimpl PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
set_target_properties(modelimpl PROPERTIES PREFIX "")
