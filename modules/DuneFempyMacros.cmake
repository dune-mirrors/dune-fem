find_package(PythonInterp 2.7 REQUIRED)
find_package(PythonLibs 2.7 REQUIRED)

include_directories("${PYTHON_INCLUDE_DIRS}")

function(add_python_target file)
  set(OUT ${CMAKE_CURRENT_BINARY_DIR}/${file}.pyc)
  set(IN ${CMAKE_CURRENT_SOURCE_DIR}/${file}.py)
  add_custom_command(OUTPUT ${OUT}
    DEPENDS ${IN}
    COMMAND ${PYTHON_EXECUTABLE} -c "import py_compile; py_compile.compile('${IN}',cfile='${OUT}')"
    VERBATIM
  )
  add_custom_target(${file}.pyc ALL DEPENDS ${OUT})
endfunction()
