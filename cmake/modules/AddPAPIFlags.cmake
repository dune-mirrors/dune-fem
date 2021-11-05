# Module providing convenience functions for using 
#
# Provides the following functions:
#
# add_dune_papi_flags(target1 target2 ...)
#
# Adds the necessary flags to compile and link the targets with PAPI support.
#
function(add_dune_papi_flags _targets)
  if(PAPI_FOUND)
    
    foreach(_target ${_targets})
      target_link_libraries(${_target} PUBLIC ${PAPI_LIBRARIES})
    endforeach(_target ${_targets})
    set_property(TARGET ${_target}
      APPEND_STRING
      PROPERTY COMPILE_FLAGS ENABLE_PAPI=1 )
    set_property(TARGET ${_targets} APPEND PROPERTY
      INCLUDE_DIRECTORIES "${PAPI_INCLUDE_DIRS}")
  endif(PAPI_FOUND)
endfunction(add_dune_papi_flags _targets)
