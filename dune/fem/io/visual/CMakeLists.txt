set(AUTOMAKE_OPTIONS  foreign 1.5)
set(SUBDIRS  grape)
set(visualdir   ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/io/visual)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt EXCLUDE_FROM_ALL)
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
endforeach(i ${SUBDIRS})
