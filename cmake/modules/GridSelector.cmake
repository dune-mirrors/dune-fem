function( get_dune_gridselector_target )
  string( REPLACE "${CMAKE_BINARY_DIR}" "" relative_dir "${CMAKE_CURRENT_BINARY_DIR}" )
  string( REPLACE "/" "_" relative_target ${relative_dir} )
  set( dune_grid_selector_target "dune_grid_selector${relative_target}" PARENT_SCOPE )
endfunction()

function( add_dune_gridselector Output )
  if( ${ARGC} GREATER 1 )
    message( SEND_ERROR "add_dune_gridselector( <Output> )" )
  endif()

  if( ${CMAKE_PROJECT_NAME} STREQUAL "dune-fem" )
    set( DUNE_GRID_SELECT_COMMAND ${CMAKE_SOURCE_DIR}/bin/dune-grid-select )
  else()
    set( DUNE_GRID_SELECT_COMMAND ${dune-fem_PREFIX}/bin/dune-grid-select )
  endif()

  set( DUNE_GRID_SELECTOR_INCLUDES )
  foreach( definition ${DUNE_GRID_SELECTOR_DEFINITIONS} )
    set( DUNE_GRID_SELECTOR_INCLUDES "${DUNE_GRID_SELECTOR_INCLUDES} -I ${definition}" )
  endforeach()

  if( UNIX )
    file( WRITE ${CMAKE_CURRENT_BINARY_DIR}/dune-grid-select
      "${DUNE_GRID_SELECT_COMMAND}${DUNE_GRID_SELECTOR_INCLUDES} -o ${CMAKE_CURRENT_BINARY_DIR}/${Output} \"\$@\""
    )
    execute_process( COMMAND chmod "+x" "${CMAKE_CURRENT_BINARY_DIR}/dune-grid-select" )
  endif()

  string( REGEX REPLACE "[.]hh$" ".dgs" Input ${Output} )
  set( DUNE_GRID_SELECTOR_INCLUDES )
  foreach( definition ${DUNE_GRID_SELECTOR_DEFINITIONS} )
    set( DUNE_GRID_SELECTOR_INCLUDES ${DUNE_GRID_SELECTOR_INCLUDES} -I ${definition} )
  endforeach()
  add_custom_command(
    OUTPUT ${Output}
    DEPENDS ${Input}
    COMMAND ${DUNE_GRID_SELECT_COMMAND}
    ARGS ${DUNE_GRID_SELECTOR_INCLUDES} -o ${CMAKE_CURRENT_BINARY_DIR}/${Output} -f ${Input}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    VERBATIM
  )

  get_dune_gridselector_target()
  add_custom_target( ${dune_grid_selector_target} DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${Output} )
endfunction()

function( target_depend_dune_gridselector Target )
  get_dune_gridselector_target()
  add_dependencies( ${Target} ${dune_grid_selector_target} )
  target_include_directories( ${Target} PRIVATE ${CMAKE_CURRENT_BINARY_DIR} )
endfunction()

macro( add_dune_grid_selector_definitions Module )
  if( ${ARGC} EQUAL 1 )
    set( ADD_DUNE_GRID_SELECTOR_DEFINITIONS ${Module} )
  else()
    set( ADD_DUNE_GRID_SELECTOR_DEFINITIONS ${ARGN} )
  endif()
  if( ${Module} STREQUAL ${CMAKE_PROJECT_NAME} )
    foreach( Definition ${ADD_DUNE_GRID_SELECTOR_DEFINITIONS} )
      list( APPEND DUNE_GRID_SELECTOR_DEFINITIONS ${CMAKE_SOURCE_DIR}/gridselector/${Definition} )
    endforeach()
  else()
    foreach( Definition ${ADD_DUNE_GRID_SELECTOR_DEFINITIONS} )
      list( APPEND DUNE_GRID_SELECTOR_DEFINITIONS ${${Module}_PREFIX}/gridselector/${Definition} )
    endforeach()
  endif()
endmacro()
