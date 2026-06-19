# FindPETSc4py.cmake
execute_process(COMMAND ${DUNE_PYTHON_VIRTUALENV_EXECUTABLE} -m pip list)
execute_process(COMMAND ${DUNE_PYTHON_VIRTUALENV_EXECUTABLE} -c
  "import petsc4py,sys; print(petsc4py.get_include())"
  OUTPUT_VARIABLE PETSC4PY_INCLUDE OUTPUT_STRIP_TRAILING_WHITESPACE)
if(PETSC4PY_INCLUDE)
  set(PETSC4PY_INCLUDES ${PETSC4PY_INCLUDE} CACHE PATH "petsc4py include")
  mark_as_advanced(PETSC4PY_INCLUDES)
  set(PETSC4PY_FOUND TRUE)
else()
  set(PETSC4PY_FOUND FALSE)
file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
  "Determing location of PETSc4py failed:\n"
  "using Python ${DUNE_PYTHON_VIRTUALENV_EXECUTABLE}\n\n")
endif()
