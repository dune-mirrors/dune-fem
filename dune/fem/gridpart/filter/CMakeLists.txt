set(filterdir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/gridpart/filter)
set(SUBDIRS )
set(filter_HEADERS  basicfilterwrapper.hh  filter.hh
                 inversefilter.hh  radialfilter.hh
                 threadfilter.hh  voidfilter.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${filter_HEADERS} DESTINATION ${filterdir})
foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt EXCLUDE_FROM_ALL)
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
endforeach(i ${SUBDIRS})
