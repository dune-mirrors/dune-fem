# some functions supporting some abbreviations

# usage:
# 
# dune_add_subdirs( <dir1> [<dir2>...]
#                   [EXCLUDE <exclude string> | NOEXCLUDE]
#
# adds for each <dir1>, <dir2>... a subdirectory except the subdirectory 
# name is equal to <exclude string>.
#
# If no <exclude string> is given a default value "test" will be used.
#
# Neglect the default value by using the NOEXCLUDE option.
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

# usage
#
# dune_add_test( <test1> [<test2>...] 
#                [TEST <test_1> [test_2] ]
#                [DEPENDENCY <dep1> [dep_2] ]
#
# Simply adds all test to the cmake testing framework.
#
# All of these tests are not build during make all.
# Build these targets with make test.
function(dune_add_test)
  include(CMakeParseArguments)
  cmake_parse_arguments(ADD_TEST "" "" "TEST;DEPENDENCY" ${ARGN})
    
  foreach(i ${ADD_TEST_UNPARSED_ARGUMENTS} ${ADD_TEST_TEST})
    add_test(${i} ${i})
  endforeach(i ${ADD_TEST_UNPARSED_ARGUMENTS} ${ADD_TEST_TEST})

  if(ADD_TEST_UNPARSED_ARGUMENTS)
    add_directory_test_target(_test_target)
    add_dependencies(${_test_target} ${ADD_TEST_UNPARSED_ARGUMENTS} ${ADD_TEST_DEPENDENCY} ${COMPILE_XFAIL_TESTS})
  endif()
endfunction(dune_add_test)
