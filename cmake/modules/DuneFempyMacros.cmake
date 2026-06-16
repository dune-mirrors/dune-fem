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


find_package(PETSc4py)

if(PETSC4PY_FOUND)
  # set HAVE_PETSC4PY for config.h
  set(HAVE_PETSC4PY ${PETSC4PY_FOUND})
  # log result
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
    "Determining location of PETSC4PY succeded:\n"
    "Include directory: ${PETSC4PY_INCLUDES}\n\n")
  set(PETSC4PY_DUNE_COMPILE_FLAGS "-DENABLE_PETSC4PY=1")
  foreach(dir ${PETSC4PY_INCLUDES})
    set(PETSC4PY_DUNE_COMPILE_FLAGS "${PETSC4PY_DUNE_COMPILE_FLAGS} -I${dir}")
  endforeach()
  set(PETSC4PY_DUNE_COMPILE_FLAGS "${PETSC4PY_DUNE_COMPILE_FLAGS}"
    CACHE STRING "Compile Flags used by DUNE when compiling with PETSc4py programs")
  foreach(dir ${PETSC4PY_INCLUDES})
    set_property(GLOBAL APPEND PROPERTY ALL_PKG_FLAGS "-I${dir}")
  endforeach()
  # register includes and libs to global flags
  dune_register_package_flags(COMPILE_DEFINITIONS "ENABLE_PETSC4PY=1"
                              INCLUDE_DIRS ${PETSC4PY_INCLUDES})
else()
  # log errornous result
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
    "Determing location of PETSc4py failed:\n"
    "Include directory: ${PETSC4PY_INCLUDES}\n\n")
endif()
