# set up project
project("dune-fem" C CXX)

#circumvent not building docs
set(BUILD_DOCS 1)
set(DUNE_REENABLE_ADD_TEST TRUE)

# general stuff
cmake_minimum_required(VERSION 2.8.12)

if(NOT (dune-common_DIR OR dune-common_ROOT OR
        "${CMAKE_PREFIX_PATH}" MATCHES ".*dune-common.*"))
    string(REPLACE  ${CMAKE_PROJECT_NAME} dune-common dune-common_DIR
      ${PROJECT_BINARY_DIR})
endif()

#find dune-common and set the module path
find_package(dune-common REQUIRED)
list(APPEND CMAKE_MODULE_PATH ${dune-common_MODULE_PATH})

#include the dune macros
include(DuneMacros)

# start a dune project with information from dune.module
dune_project()
dune_enable_all_packages()

include(DuneMPI)

#add sub_directories
dune_add_subdirs(dune lib doc cmake/modules)

add_subdirectory("python")
dune_python_install_package(PATH "python")

# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
