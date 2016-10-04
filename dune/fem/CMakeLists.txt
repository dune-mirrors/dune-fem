add_subdirectory(space)
add_subdirectory(discretefunction)
add_subdirectory(scheme)
add_subdirectory(function)
add_subdirectory(model)

add_library(_fem SHARED _fem.cc)
set_target_properties(_fem PROPERTIES PREFIX "")

add_python_targets(fem
  __init__
  _adaptation
  view
)
