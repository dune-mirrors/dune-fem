# some functions supporting some command line hacks
#



# makes it possible to define GRIDDIM, WORLDDIM or GRIDTYPE
# via command line (or export function)
# (example: make GRIDDIM=2 WORLDDIM=2 GRIDTYPE=SGRID _targets)
#
#
# Call the set_gridtype_from_command_line()-function to use this command line option
# in the following way
#
# set_gridtype_from_command_line(<target> [<target2>...] 
#                                [GRIDTYPE <gridtype>]
#                                [GRIDDIM <griddim>]
#                                [WORLDDIM <worlddim>]  )
#
# or use the second version:
#
# set_gridtype_from_command_line(<target> [<target2>...] 
#                                 [FORCE] )
#
# GRIDTYPE <arg>:     use the argument <arg> instead of command line option
# GRIDDIM <arg>:      use the argument <arg> instead of command line option
# WORLDDIM <arg>:     use the argument <arg> instead of command line option
#
# FORCE:        enable force option to force reading all three arguments via command line:
#               You can neglect the other parameters GRIDTYPE, GRIDDIM and WOLRDDIM
#
# Notice: Once, GRIDDIM, WORLDDIM or GRIDTYPE can be chosen via command line
# these compile definitions _HAVE_ to be specified via command line. 
#
# Warning: Use this function with caution! It may be broken on some systems
function(set_gridtype_from_command_line)

  include(CMakeParseArguments)
  cmake_parse_arguments(CMDLINE_PARAM "FORCE" "GRIDTYPE;GRIDDIM;WORLDDIM" "" ${ARGN})

  #known grids
  if( NOT GRID_LIST)
    create_grid_list()
  endif()

  # this is a real dirty hack. 
  # At least "=$(GRIDTYPE)" should fail on some cmake generators/systems, I guess...
  foreach(_targets ${CMDLINE_PARAM_UNPARSED_ARGUMENTS})

    #erase existing 
    get_target_property(_props_compile_definitions ${_targets} COMPILE_DEFINITIONS)
    string(REPLACE "_props_compile_definitions-NOTFOUND" "" _props_compile_definitions "${_props_compile_definitions}")

    #delete some of existing flags 
    set(compile_def ${_props_compile_definitions})
    if( compile_def )
      foreach(_del_grids ${GRID_LIST})
        list(REMOVE_ITEM compile_def "${_del_grids}")
      endforeach(_del_grids ${GRID_LIST})
    endif()
    set(compile_def_save ${compile_def})
    foreach(_flags_ ${compile_def_save})
      string(REGEX MATCH "GRIDDIM=[0123456789]?" _result ${_flags_} )
      list(REMOVE_ITEM compile_def "${_result}")
      string(REGEX MATCH "WORLDDIM=[0123456789]?" _result ${_flags_} )
      list(REMOVE_ITEM compile_def "${_result}")
    endforeach()
    set_property(TARGET ${_targets} PROPERTY COMPILE_DEFINITIONS "${compile_def}")

    if( CMDLINE_PARAM_FORCE )
      set_property(TARGET ${_targets} APPEND PROPERTY 
        COMPILE_DEFINITIONS "WORLDDIM=$(WORLDDIM);GRIDDIM=$(GRIDDIM);=$(GRIDTYPE)" )
    else()
      if( CMDLINE_PARAM_GRIDDIM )
        set_property(TARGET ${_targets} APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=${CMDLINE_PARAM_GRIDDIM}" )
        if( CMDLINE_PARAM_WORLDDIM )
          set_property(TARGET ${_targets} APPEND PROPERTY COMPILE_DEFINITIONS "WORLDDIM=${CMDLINE_PARAM_WORLDDIM}" )
        else()
          #default value it GRIDDIM is given
          set_property(TARGET ${_targets} APPEND PROPERTY COMPILE_DEFINITIONS "WORLDDIM=${CMDLINE_PARAM_GRIDDIM}" )
        endif()
      else()
        set_property(TARGET ${_targets} APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=$(GRIDDIM)" )
        if( CMDLINE_PARAM_WORLDDIM )
          set_property(TARGET ${_targets} APPEND PROPERTY COMPILE_DEFINITIONS "WORLDDIM=${CMDLINE_PARAM_WORLDDIM}" )
        else()
          set_property(TARGET ${_targets} APPEND PROPERTY COMPILE_DEFINITIONS "WORLDDIM=$(WORLDDIM)" )
        endif()
      endif()
      if( CMDLINE_PARAM_GRIDTYPE )
        set_property(TARGET ${_targets} APPEND PROPERTY COMPILE_DEFINITIONS "${CMDLINE_PARAM_GRIDTYPE}" )
      else()
        set_property(TARGET ${_targets} APPEND PROPERTY COMPILE_DEFINITIONS "=$(GRIDTYPE)" )
      endif()
    endif()
  endforeach()  

endfunction(set_gridtype_from_command_line _targets)



function(create_grid_list)
  #known grids
  set( GRID_LIST SGRID YASPGRID ONEDGRID CACHE STRING 
       "grid compile definitions known in dune")
  if( HAVE_ALUGRID ) 
    set( GRID_LIST ${GRID_LIST} ALUGRID_CUBE ALUGRID_CONFORM ALUGRID_SIMPLEX )
  endif()
  if( HAVE_ALBERTA )
    set( GRID_LIST ${GRID_LIST} ALBERTAGRID )
  endif()
  if( HAVE_UG ) 
    set( GRID_LIST ${GRID_LIST} UGGRID )
  endif()  
endfunction(create_grid_list)

# adds for all/many GRIDTYPE/GRIDDIM/WORLDDIM combinations executables
# 
# Use this function in the following way:
#
# dune_create_gridtype_executable(<target> [target2>...]
#                                 [WIN32][MACOSX_BUNDLE][EXCLUDE_FROM_ALL] 
#                                 )
# 
# WIN32 and MACOSX_BUNDLE are the same as for add_executable().
# By default all generated are targets are "EXCLUDE_FROM_ALL".
#
# Notice: We suppose that all flags etc. are set up correctly for the 
# given targets <target>, <target2>...
# 
# The new executables are name in the following way:
# <target>_<gridtype>_<griddim>_<worlddim>
function(dune_create_gridtype_executable)

  
  include(CMakeParseArguments)
  cmake_parse_arguments(ADD_EXEC "WIN32;MACOSX_BUNDLE;EXCLUDE_FROM_ALL" "" "" ${ARGN}) 
  
  foreach(_targets ${ARGN})

  
    #we expect, that the given executable is set up correctly...
    get_target_property(_props_compile_flags ${_targets} COMPILE_FLAG)
    get_target_property(_props_compile_definitions ${_targets} COMPILE_DEFINITIONS)
    get_target_property(_props_include_directories ${_targets} INCLUDE_DIRECTORIES)
    get_target_property(_props_link_libraries ${_targets} LINK_LIBRARIES)
    get_target_property(_props_sources ${_targets} SOURCES)
  
    string(REPLACE "_props_compile_flags-NOTFOUND" "" _props_compile_flags "${_props_compile_flags}")
    string(REPLACE "_props_compile_definitions-NOTFOUND" "" _props_compile_definitions "${_props_compile_definitions}")
    string(REPLACE "_props_include_directories-NOTFOUND" "" _props_include_directories "${_props_include_directories}")
    string(REPLACE "_props_link_libraries-NOTFOUND" "" _props_link_libararies "${_props_link_libraries}")
    string(REPLACE "_props_sources-NOTFOUND" "" _props_sources "${_props_sources}")
  
    #known grids
    if( NOT GRID_LIST)
      create_grid_list()
    endif()

    #known grids with dimensions... 
    set( GRID_LIST_ALL "ONEDGRID 1 1"
                       "YASPGRID 1 1"
                       "YASPGRID 2 2"
                       "YASPGRID 3 3"
                       "SGRID 1 1"
                       "SGRID 2 2"
                       "SGRID 3 3" 
                       )
    if( HAVE_ALUGRID ) 
      set( GRID_LIST_ALL ${GRID_LIST_ALL}
                         "ALUGRID_CUBE 1 1"
                         "ALUGRID_CUBE 2 2"
                         "ALUGRID_CUBE 3 3"
                         "ALUGRID_SIMPLEX 1 1"
                         "ALUGRID_SIMPLEX 2 2"
                         "ALUGRID_SIMPLEX 3 3"
                         "ALUGRID_CONFORM 1 1"
                         "ALUGRID_CONFORM 2 2"
                         "ALUGRID_CONFORM 3 3" 
                         )
    endif()
    if( HAVE_ALBERTA )
      set( GRID_LIST_ALL ${GRID_LIST_ALL}
                         "ALBERTA 1 1"
                         "ALBERTA 2 2"
                         "ALBERTA 3 3"
                         )
    endif()
    if( HAVE_UG ) 
      set( GRID_LIST_ALL ${GRID_LIST_ALL}
                         "UGGRID 1 1"
                         "UGGRID 2 2"
                         "UGGRID 3 3"
                         )
    endif()
    
    #delete some of existing flags 
    set(compile_flags ${_props_compile_definitions})
    if( compile_flags )
      foreach(_del_grids ${GRID_LIST})
        list(REMOVE_ITEM compile_flags "${_del_grids}")
      endforeach(_del_grids ${GRID_LIST})
    endif()
    set(compile_flags_save ${compile_flags})
    foreach(_flags_ ${compile_flags_save})
      string(REGEX MATCH "GRIDDIM=[0123456789]?" _result ${_flags_} )
      list(REMOVE_ITEM compile_flags "${_result}")
      string(REGEX MATCH "WORLDDIM=[0123456789]?" _result ${_flags_} )
      list(REMOVE_ITEM compile_flags "${_result}")
    endforeach()
  
    #build all targets
    foreach(_grids ${GRID_LIST_ALL} )
      separate_arguments(_grids_sep UNIX_COMMAND "${_grids}")
      list(GET _grids_sep 0 _gridtype)
      list(GET _grids_sep 1 _griddim)
      list(GET _grids_sep 2 _worlddim)
   
      #new target name 
      set(_target_name ${_targets}_${_gridtype}_${_griddim}_${_worlddim} )
  
      #append new compile flags regarding gridtype, grid dim and world dim
      set(cur_compile_flags "GRIDDIM=${_griddim};WORLDDIM=${_worlddim};${_gridtype};${compile_flags}")
      
      #create executable and add all important flags...
      if( ADD_EXEC_MACOSX_BUNDLE )
        add_executable(${_target_name} MACOSX_BUNDLE EXCLUDE_FROM_ALL ${_props_sources}) 
      elseif(ADD_EXEC_WIN32)
        add_executable(${_target_name} WIN32 EXCLUDE_FROM_ALL ${_props_sources}) 
      else()
        add_executable(${_target_name} EXCLUDE_FROM_ALL ${_props_sources}) 
      endif()

      set_property(TARGET ${_target_name} APPEND PROPERTY
        LINK_LIBRARIES "${_props_link_libraries}")
      set_property(TARGET ${_target_name} APPEND PROPERTY 
        COMPILE_DEFINITIONS "${compile_definitions}" )
      set_property(TARGET ${_target_name} APPEND PROPERTY
        INCLUDE_DIRECTORIES "${_props_include_directories}")
      set_property(TARGET ${_target_name} APPEND PROPERTY
        COMPILE_FLAGS "${cur_compile_flags}")
    endforeach()
  
  endforeach()
endfunction(dune_create_gridtype_executable)
