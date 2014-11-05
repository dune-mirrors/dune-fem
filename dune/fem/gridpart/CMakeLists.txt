set(gridpartdir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/gridpart)
set(SUBDIRS  common  filter  filteredgridpart  geogridpart  idgridpart  test)
set(gridpart_HEADERS adaptiveleafgridpart.hh
                     adaptiveleafindexset.hh
                     codimindexset.hh
                     defaultindexsets.hh
                     dunefemindexsets.hh
                     emptyindexset.hh
                     filteredgridpart.hh
                     filter.hh
                     geogridpart.hh
                     idgridpart.hh
                     leafgridpart.hh
                     levelgridpart.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${gridpart_HEADERS} DESTINATION ${gridpartdir})
dune_add_subdirs(${SUBDIRS})
