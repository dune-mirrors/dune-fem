add_subdirectory(gridpart)
add_subdirectory(space)
add_subdirectory(discretefunction)
add_subdirectory(scheme)

add_library(femmpi SHARED femmpi.cc)
set_target_properties(femmpi PROPERTIES PREFIX "")

# add_python_target(grid)
# add_python_target(create)
# add_python_target(function)

# configure_file(__init__.py __init__.py COPYONLY)
