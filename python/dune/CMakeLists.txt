add_subdirectory(fem)
add_subdirectory(fenics)
add_subdirectory(models)
add_subdirectory(source)
add_subdirectory(ufl)

add_python_targets(dune
  __init__
)
