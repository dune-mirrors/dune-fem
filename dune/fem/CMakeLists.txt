add_subdirectory(gridpart)
add_subdirectory(space)
add_subdirectory(discretefunction)
add_subdirectory(scheme)

add_library(femmpi SHARED femmpi.cc)
set_target_properties(femmpi PROPERTIES PREFIX "")

add_python_targets(fem
  __init__
  grid
  create
  function
)
