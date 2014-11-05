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
dune_add_subdirs(${SUBDIRS})

dune_add_library(io OBJECT io.cc)
