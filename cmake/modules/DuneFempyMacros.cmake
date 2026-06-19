# find python package ufl with at least version 2016.2
# required cannot be used otherwise the version is fixed to 2016.2
dune_python_find_package(PACKAGE "ufl" VERSION "2016.2")

include(UseModelCompiler)

function(add_gmshgeo_target file)
  set(OUT ${CMAKE_CURRENT_BINARY_DIR}/${file}.msh)
  set(IN ${CMAKE_CURRENT_SOURCE_DIR}/${file}.geo)
  add_custom_command(OUTPUT ${OUT}
    DEPENDS ${IN}
    COMMAND gmsh -2 -o ${OUT} ${IN}
    VERBATIM
  )
  add_custom_target(${file}.msh ALL DEPENDS ${OUT})
endfunction()
