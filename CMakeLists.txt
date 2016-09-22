add_subdirectory(dune)

find_program(PYTHON "python")
if (PYTHON)
    file(GLOB_RECURSE MY_PYS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "dune/*.py")
    foreach(ITEM ${MY_PYS})
      # get_filename_component(ITEM_PATH ${ITEM} PATH)
      # file(COPY ${ITEM} DESTINATION ${ITEM_PATH})
      configure_file(${ITEM} ${ITEM} COPYONLY)
    endforeach()

    set(SETUP_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in")
    set(SETUP_PY    "${CMAKE_CURRENT_BINARY_DIR}/setup.py")
    configure_file(${SETUP_PY_IN} ${SETUP_PY})

    if(PYTHON_INSTALL_LOCATION STREQUAL "user")
      set(PYTHON_INSTALL_LOCATION "--user")
      message("--user")
    elseif(PYTHON_INSTALL_LOCATION STREQUAL "global")
      set(PYTHON_INSTALL_LOCATION "")
      message("--global")
    else()
      set(PYTHON_INSTALL_LOCATION "--prefix=${CMAKE_INSTALL_PREFIX}")
      message("--cmake")
    endif()
    install(CODE
      "execute_process(COMMAND python setup.py install -f ${PYTHON_INSTALL_LOCATION}
               WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})")
endif()
