add_subdirectory(discretefunction)
add_subdirectory(function)
add_subdirectory(model)
add_subdirectory(operator)
add_subdirectory(scheme)
add_subdirectory(space)
add_subdirectory(utility)

add_python_targets(fem
  __init__
  __main__
  _adaptation
  _spaceadaptation
  plotting
  view
)

dune_add_pybind11_module(NAME _fem)
set_property(TARGET _fem PROPERTY LINK_LIBRARIES dunecommon dunegeometry dunegrid dunefem APPEND)
if(SKBUILD)
  install(TARGETS _fem LIBRARY DESTINATION python/dune/fem)
endif()
