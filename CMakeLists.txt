# set up project
project("dune-fem" C CXX)

#circumvent not building docs
set(BUILD_DOCS 1)

# general stuff
cmake_minimum_required(VERSION 2.8.6)

#find dune-common and set the module path
find_package(dune-common)
list(APPEND CMAKE_MODULE_PATH ${dune-common_MODULE_PATH}
     "${CMAKE_SOURCE_DIR}/cmake/modules")# make sure our own modules are found

#include the dune macros
include(DuneMacros)

# start a dune project with information from dune.module
dune_project()

set(SUBDIRS  dune lib m4 doc)
# we need the module file to be able to build via dunecontrol
set(EXTRA_DIST  dune.module ChangeLog modules)
set(DISTCHECK_CONFIGURE_FLAGS  --with-dune="${DUNEROOT}" CXX="${CXX}" CC="${CC}")
foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt EXCLUDE_FROM_ALL)
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
endforeach(i ${SUBDIRS})

add_subdirectory(cmake/modules)
# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
