# Module providing convenience functions for using
#
# Provides the following functions:
#
# add_dune_psctoolkit_flags(target1 target2 ...)
#
# Adds the necessary flags to compile and link the targets with PSCTOOLKIT support.
#
function(add_dune_psctoolkit_flags _targets)
  if(PSCTOOLKIT_FOUND)
    foreach(_target ${_targets})
      target_link_libraries(${_target} PUBLIC ${PSCTOOLKIT_LIBRARIES})
    endforeach(_target ${_targets})
    set_property(TARGET ${_target}
      APPEND_STRING
      PROPERTY COMPILE_FLAGS ENABLE_PSCTOOLKIT=1 )
    set_property(TARGET ${_targets} APPEND PROPERTY
      INCLUDE_DIRECTORIES "${PSCTOOLKIT_INCLUDE_DIRS}")
  endif(PSCTOOLKIT_FOUND)
endfunction(add_dune_psctoolkit_flags _targets)
