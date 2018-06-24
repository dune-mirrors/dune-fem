set(DUNE_BUILD_SET_PYTHON_PATH "import sys")
foreach(module ${ALL_DEPENDENCIES})
  if(EXISTS "${${module}_DIR}/python" AND IS_DIRECTORY "${${module}_DIR}/python")
    list(APPEND DUNE_BUILD_SET_PYTHON_PATH "sys.path.append('${${module}_DIR}/python')")
  endif()
endforeach()

if(EXISTS "${CMAKE_SOURCE_DIR}/python" AND IS_DIRECTORY "${CMAKE_SOURCE_DIR}/python")
  list(APPEND DUNE_BUILD_SET_PYTHON_PATH "sys.path.append('${CMAKE_BINARY_DIR}/python')")
endif()

function(dune_fem_add_elliptic_models)
  foreach(source ${ARGN})
    get_filename_component(base ${source} NAME_WE)
    add_custom_command(
        OUTPUT ${base}.hh
        DEPENDS ${source}
        COMMAND ${PYTHON_EXECUTABLE} -c "${DUNE_BUILD_SET_PYTHON_PATH}; from dune.models.elliptic.compiler import main; sys.exit(main())" ${source} -o ${CMAKE_CURRENT_BINARY_DIR}/${base}.hh -n ${base}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMENT "Compiling elliptic model ${base}..."
        VERBATIM
      )
    set_property(SOURCE ${CMAKE_CURRENT_BINARY_DIR}/${base}.hh PROPERTY GENERATED 1)
  endforeach()
endfunction()

function(dune_fem_add_integrands_models)
  foreach(source ${ARGN})
    get_filename_component(base ${source} NAME_WE)
    add_custom_command(
        OUTPUT ${base}.hh
        DEPENDS ${source}
        COMMAND ${PYTHON_EXECUTABLE} -c "${DUNE_BUILD_SET_PYTHON_PATH}; from dune.models.integrands.compiler import main; sys.exit(main())" ${source} -o ${CMAKE_CURRENT_BINARY_DIR}/${base}.hh
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMENT "Compiling integrands model ${base}..."
        VERBATIM
      )
    set_property(SOURCE ${CMAKE_CURRENT_BINARY_DIR}/${base}.hh PROPERTY GENERATED 1)
  endforeach()
endfunction()
