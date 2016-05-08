add_subdirectory(fem)
add_subdirectory(generated)
add_subdirectory(generator)

add_library(common SHARED common.cc)
set_target_properties(common PROPERTIES PREFIX "")

add_library(femmpi SHARED femmpi.cc)
set_target_properties(femmpi PROPERTIES PREFIX "")

configure_file(__init__.py __init__.py COPYONLY)
