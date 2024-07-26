# general stuff
cmake_minimum_required(VERSION 3.16)

# set up project
project("dune-fem" C CXX)

#circumvent not building docs
set(BUILD_DOCS 1)
set(DUNE_REENABLE_ADD_TEST TRUE)

# if DUNE_REQUIRE_PYTHONBINDINGS was NOT explicitly disabled then Python is required!
if(DEFINED DUNE_REQUIRE_PYTHONBINDINGS AND NOT DUNE_REQUIRE_PYTHONBINDINGS)
  set(DUNE_PYTHONBINDINGS_REQUIRED FALSE)
else()
  set(DUNE_PYTHONBINDINGS_REQUIRED TRUE)
endif()

if(NOT (dune-common_DIR OR dune-common_ROOT OR
        "${CMAKE_PREFIX_PATH}" MATCHES ".*dune-common.*"))
    string(REPLACE  ${CMAKE_PROJECT_NAME} dune-common dune-common_DIR
      ${PROJECT_BINARY_DIR})
endif()

#find dune-common and set the module path
find_package(dune-common REQUIRED)
list(APPEND CMAKE_MODULE_PATH ${dune-common_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake/modules)

#include the dune macros
include(DuneMacros)

# start a dune project with information from dune.module
dune_project()
dune_enable_all_packages()

# Create the dune-fem library
dune_add_library(dunefem
  EXPORT_NAME Fem
  LINK_LIBRARIES ${DUNE_LIBS})

include(AddMPIFlags)

# if Python bindings are enabled, include necessary sub directories.
if( DUNE_ENABLE_PYTHONBINDINGS )
  set(PYDEMO_DIR pydemo python)
endif()

#add sub_directories
dune_add_subdirs(dune lib doc cmake/modules ${PYDEMO_DIR})

# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project()
