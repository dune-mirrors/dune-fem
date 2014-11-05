# some functions supporting some abbreviations
#

function(dune_add_subdirs)

  include(CMakeParseArguments)
  cmake_parse_arguments(SHORT "NOEXCLUDE" "EXCLUDE" "" ${ARGN})

  if( SHORT_NOEXCLUDE )
    foreach(i ${SHORT_UNPARSED_ARGUMENTS})
      add_subdirectory(${i})
    endforeach(i ${SHORT_UNPARSED_ARGUMENTS})
  else()
    if(NOT SHORT_EXCLUDE)
      #default
      set(SHORT_EXCLUDE "test")
    endif()
    foreach(i ${SHORT_UNPARSED_ARGUMENTS})
      if(${i} STREQUAL "${SHORT_EXCLUDE}")
        set(opt EXCLUDE_FROM_ALL)
      endif(${i} STREQUAL "${SHORT_EXCLUDE}")
      add_subdirectory(${i} ${opt})
    endforeach(i ${SHORT_UNPARSED_ARGUMENTS})
  endif()


endfunction(dune_add_subdirs)

