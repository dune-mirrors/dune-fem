add_subdirectory(fem)
add_subdirectory(generated)
add_subdirectory(models)

add_library(femmpi SHARED femmpi.cc)
set_target_properties(femmpi PROPERTIES PREFIX "")

configure_file(dune-fempy-setup.py dune-fempy-setup.py COPYONLY)

add_python_target(source)
add_python_target(ufl)
