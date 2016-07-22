add_subdirectory(fem)
add_subdirectory(generated)
add_subdirectory(models)

# add_library(common SHARED common.cc)
# set_target_properties(common PROPERTIES PREFIX "")

add_library(femmpi SHARED femmpi.cc)
set_target_properties(femmpi PROPERTIES PREFIX "")

add_python_target(source)
add_python_target(ufl)
configure_file(__init__.py.in __init__.py)
