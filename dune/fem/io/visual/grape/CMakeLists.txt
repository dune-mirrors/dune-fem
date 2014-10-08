set(SUBDIRS  datadisp)
set(grapedir   ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/io/visual/grape)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt EXCLUDE_FROM_ALL)
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
endforeach(i ${SUBDIRS})
