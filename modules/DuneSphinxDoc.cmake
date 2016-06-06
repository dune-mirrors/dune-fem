find_package(Sphinx)
find_package(PythonInterp)

function(dune_sphinx_doc)
  # Only proceed if Sphinx was found on the system
  if(NOT SPHINX_FOUND)
    message("-- Skipping building Sphinx documentation (Sphinx was not found!)")
    return()
  endif()

  # Only proceed if the python interpreter was found by cmake
  if(NOT PYTHONINTERP_FOUND)
    message("-- Skipping building Sphinx documentation (Python interpreter was not found!)")
    return()
  endif()

  # Parse Arguments
  include(CMakeParseArguments)
  cmake_parse_arguments(SPHINX_DOC "" "CONF" "BUILDTYPE" ${ARGN})
  if(SPHINX_DOC_UNPARSED_ARGUMENTS)
    message(WARNING "Unparsed arguments in dune_sphinx_doc")
  endif()

  # copy conf.py into build directory
  if(NOT SPHINX_DOC_CONF)
    set(SPHINX_DOC_CONF conf.py)
  endif()
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${SPHINX_DOC_CONF}.in)
    configure_file(${SPHINX_DOC_CONF}.in ${CMAKE_CURRENT_BINARY_DIR}/conf.py)
  elseif(EXISTS ${CMAKE_CURRENT_SOUREC_DIR}/${SPHINX_DOC_CONF})
    configure_file(${SPHINX_DOC_CONF} ${CMAKE_CURRENT_BINARY_DIR}/conf.py COPYONLY)
  else()
    message(SEND_ERROR "Sphinx configuration '${SPHINX_DOC_CONF}' not found.")
  endif()

  # call Sphinx for each requested build type
  if(NOT SPHINX_DOC_BUILDTYPE)
    set(SPHINX_DOC_BUILDTYPE html)
  endif()
  foreach(type ${SPHINX_DOC_BUILDTYPE})
    add_custom_target(sphinx_doc_${type}
                      COMMAND ${SPHINX_EXECUTABLE}
                                -b ${type}
                                -w ${CMAKE_BINARY_DIR}/Sphinx-${type}.log
                                -c ${CMAKE_CURRENT_BINARY_DIR}
                                ${CMAKE_CURRENT_SOURCE_DIR}
                                ${CMAKE_CURRENT_BINARY_DIR}/${type}
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                      DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/conf.py
                     )
    add_dependencies(doc sphinx_doc_${type})
  endforeach()
endfunction()
