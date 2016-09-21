add_subdirectory(fem)
add_subdirectory(models)

add_library(femmpi SHARED femmpi.cc)
set_target_properties(femmpi PROPERTIES PREFIX "")

add_python_target(__init__)
add_python_target(source)
add_python_target(ufl)
