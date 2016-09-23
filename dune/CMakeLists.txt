add_subdirectory(fem)
add_subdirectory(models)

add_python_targets(dune
  __init__
  source
  ufl
)
