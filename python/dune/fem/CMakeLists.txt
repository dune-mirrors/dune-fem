add_subdirectory(discretefunction)
add_subdirectory(function)
add_subdirectory(model)
add_subdirectory(operator)
add_subdirectory(scheme)
add_subdirectory(space)
add_subdirectory(utility)


add_python_targets(fem
  __init__
  _fem
  _adaptation
  _spaceadaptation
  plotting
  view
)

dune_symlink_to_source_files(FILES _fem.cc)
