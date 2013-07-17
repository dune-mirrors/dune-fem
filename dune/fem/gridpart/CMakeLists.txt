set(gridpartdir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/gridpart)
set(SUBDIRS  common  filter  filteredgridpart  geogridpart  idgridpart  test)
set(gridpart_HEADERS  adaptiveleafgridpart.hh
                   adaptiveleafindexset.hh
                   codimindexset.hh
                   defaultindexsets.hh
                   dunefemindexsets.hh
                   emptyindexset.hh
                   filter.hh
                   filteredgridpart.hh
                   geogridpart.hh
                   gridpart.hh
                   hierarchicgridpart.hh
                   idgridpart.hh
                   leafgridpart.hh
                   levelgridpart.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${gridpart_HEADERS} DESTINATION ${gridpartdir})
foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt EXCLUDE_FROM_ALL)
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
endforeach(i ${SUBDIRS})
