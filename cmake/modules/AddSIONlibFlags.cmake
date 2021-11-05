# Module providing convenience functions for using 
#
# Provides the following functions:
#
# add_dune_sionlib_flags(target1 target2 ...)
#
# Adds the necessary flags to compile and link the targets with SIONLIB support.
#
function(add_dune_sionlib_flags _targets)
  if(SIONLIB_FOUND)
    foreach(_target ${_targets})
      target_link_libraries(${_target} PUBLIC ${SIONLIB_LIBRARIES})
    endforeach(_target ${_targets})
    set_property(TARGET ${_targets}
      APPEND_STRING
      PROPERTY COMPILE_FLAGS ENABLE_SIONLIB=1 )
    set_property(TARGET ${_targets} APPEND PROPERTY
      INCLUDE_DIRECTORIES "${SIONLIB_INCLUDE_DIRS}")
  endif(SIONLIB_FOUND)
endfunction(add_dune_sionlib_flags _targets)
