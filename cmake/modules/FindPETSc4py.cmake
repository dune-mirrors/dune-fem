# FindPETSc4py.cmake
find_package(Python COMPONENTS Interpreter REQUIRED)
execute_process(COMMAND ${Python_EXECUTABLE} -c
  "import petsc4py,sys; print(petsc4py.get_include())"
  OUTPUT_VARIABLE PETSC4PY_INCLUDE OUTPUT_STRIP_TRAILING_WHITESPACE)
if(PETSC4PY_INCLUDE)
  set(PETSC4PY_INCLUDES ${PETSC4PY_INCLUDE} CACHE PATH "petsc4py include")
  mark_as_advanced(PETSC4PY_INCLUDES)
  set(PETSC4PY_FOUND TRUE)
else()
  set(PETSC4PY_FOUND FALSE)
endif()
