# install header files
set(iodir   ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/io)
set(io_HEADERS  io.hh  parameter.hh)
# do everything in the subdirs
set(SUBDIRS  file streams visual matlab test)
# build io library

# include common rules
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${io_HEADERS} DESTINATION ${iodir})
foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt EXCLUDE_FROM_ALL)
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
endforeach(i ${SUBDIRS})

dune_add_library(io OBJECT io.cc)
